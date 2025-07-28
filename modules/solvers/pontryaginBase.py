from itertools import product
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root
from .solver import Solver

class PontryaginCompute(Solver):
    __name__ = "Метод Понтрягина"
    __res = {
        "tk":None,
        "xn":None,
        "un":None
    }
    def __init__(self, fs, bs, args=(),u_opt=None, **kwargs) -> None:
        super().__init__()
        self.eps = 1e-5
        self.num = 50
        self.method = 'hybr' 
        self.diffmethod='RK45'
        self.findBalancePoint = False
        self.phi = None
        self.isBkFun = False
        self.psi0 = args
        for key, value in kwargs.items():
            match key:
                case 'limit':
                    self.is_limit = kwargs[key]
                case 'N':
                    self.n = int(kwargs[key])
                case 'Varepsilon':
                    if value!='':
                        self.eps = float(value)
                case 'Num':
                    if value!='':
                        self.num = int(value)
                case 'Diffmethod':self.diffmethod = value
                case 'Method':self.method = value
        self.bk = bs[self.n:2*self.n]
        self.t0 = bs[-1]
        self.b0 = bs[0:self.n]
        self.fs = fs
        self.functional = fs[-1]

    def prepareData(self):
        xs = [sp.Function(f'x_{i}') for i in range(1, self.n+1)]
        psis = list(sp.symbols(f'psi1:{self.n+1}', cls=sp.Function))
        t = sp.Symbol('t')
        u = sp.Function('u')(t)
        H = sum([psi(t)*f for psi, f in zip(psis, self.fs)])-(self.functional.function)

        if self.phi:
            self.phi = sp.lambdify([x(t) for x in xs], self.phi)

        if H.diff(u, 2) == 0:
            if self.is_limit:
                if (self.functional.function.eval()) == 1:
                    u_opt = sp.sign(H.diff(u))
                else:
                    u_opt = sp.Piecewise((-1, H.diff(u)< -1), (1, H.diff(u)>1), (H.diff(u), True))
            else:
                u_opt = H.diff(u)
        else:
            if self.is_limit:
                u_opt_part = sp.solve(H.diff(u), u)[0]
                self.u_opt = sp.Piecewise((u_opt_part, sp.Abs(u_opt_part) <= 1), (1, u_opt_part > 1), (-1, u_opt_part < -1))
            else:
                self.u_opt = sp.solve(H.diff(u), u)[0]

        if H.diff(u, 2) == 0:
            if self.is_limit:
                c = sp.Symbol('c')

                self.F =  sp.lambdify((*[x(t) for x in xs], *[psi(t) for psi in self.psis]),
                     [H.subs(u, c).diff(psi(t)).subs(c, self.u_opt) for psi in self.psis] + [-H.subs(u, c).diff(x(t)).subs(c, u_opt) for x in xs],'scipy',dummify=True)
            else:
                self.F = sp.lambdify((*[x(t) for x in xs], *[psi(t) for psi in self.psis]),
                                [H.subs(u, self.u_opt).diff(psi(t)) for psi in self.psis] + [-H.subs(u, u_opt).diff(x(t)) for
                                                                                        x in xs], ('scipy', 'numpy'),
                                dummify=True)
        else:
            self.F =  sp.lambdify((*[x(t) for x in xs], *[psi(t) for psi in self.psis]),
                 [H.subs(u, u_opt).diff(psi(t)) for psi in self.psis]+[-H.subs(u, self.u_opt).diff(x(t)) for x in xs],('scipy','numpy'),dummify=True)

        self.Hf = sp.lambdify((*[x(t) for x in xs], *[psi(t) for psi in self.psis]),H.subs(u, self.u_opt))
  
    async def findMinimum(self):
        def F_(*x):
            t, x0 = x
            return self.F(*x0)
        def Nev(psi0):
            # Проверка корректности конечного времени
            if psi0[-1] < 0:
                return [1e9] * (self.n+1)  # Возвращаем большое число, если время отрицательное

            # Интегрирование системы
            tn = np.linspace(self.t0, psi0[-1], self.num)
            sol = solve_ivp(F_, [self.t0, psi0[-1]], np.concatenate((self.b0,psi0[:-1])), t_eval=tn, method=self.diffmethod, atol=self.eps)
            psiT = sol.y[:, -1]  # Последнее значение решения
            # Выделение переменных состояния
            xf = psiT[:self.n]
            # Возвращаем разность между текущим и целевым состоянием
            return np.concatenate(([i - k for i, k in zip(xf, self.bk)], [self.Hf(*self.bk, *psiT[self.n:])]))

        psi0 = self.psi0
        
        self.y0 = await root(Nev, psi0, method=self.method, tol=self.eps)['x']
        
    async def findOptimal(self):
        tn = np.linspace(self.t0, self.y0[-1], self.num)
        y_toch = await solve_ivp(self.F_, [self.t0, tn[-1]], np.concatenate((self.b0, self.y0[:-1])), t_eval=tn, method = self.diffmethod).y.T
        self.psin = y_toch[:,self.n: self.n*2]
        xn = y_toch[:,0:self.n]
        self.xk = xn[-1]

    def checkResult(self):
            assert max([abs(i-j)**2 for i,j in zip(self.xk,self.bk)])>=0.1
            
    def updateResult(self):
        self.u = sp.lambdify((*[x(self.t) for x in self.xs], *[psi(self.t) for psi in self.psis]), self.u_opt, modules='numpy')
        u = [self.u(*x, *psi) for x, psi in zip(self.xn, self.psin)]


        self.tk = self.y0[-1]
        self.__res["tk"] = self.tk
        self.__res["xn"] = self.xn
        self.__res["un"] = u
    
    @property
    def result(self):
        return self.__res