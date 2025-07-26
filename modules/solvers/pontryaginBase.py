from itertools import product
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root


class PontryaginCompute:
    def __init__(self, fs, bs, alpha1=0, alpha2=0, args=(),u_opt=None, **kwargs) -> None:
        # def is_linear(expr, vars):
        #     for xy in product(*[vars]*len(vars)):
        #             try: 
        #                 if not sp.Eq(sp.diff(expr, *xy), 0):
        #                     return False
        #             except TypeError:
        #                 return False
        #     return True
        self.eps = 1e-5
        self.num = 50
        self.method = 'hybr' 
        self.diffmethod='RK45'
        self.findBalancePoint = False
        self.phi = None
        self.isBkFun = False
        for key, value in kwargs.items():
            match key:
                case 'tkCurrent':
                    self.tkCurrent = kwargs[key]
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
                case 'Findbalancepoint':
                    if value:
                        self.findBalancePoint=True
                case 'Isbkfun':
                    if value == 'True':
                        self.isBkFun = True
                        self.phi = kwargs['phi']
                        pass
                    else:
                        self.bk = bs[self.n:2*self.n]
                        pass


        t0 = bs[-1]
        b0 = bs[0:self.n]


        xs = [sp.Function(f'x_{i}') for i in range(1, self.n+1)]
        self.psis = list(sp.symbols(f'psi1:{self.n+1}', cls=sp.Function))
        t = sp.Symbol('t')
        u = sp.Function('u')(t)
        self.functional = sp.Integral(alpha1+alpha2*sp.Pow(u, 2), (t, t0, sp.Symbol("t_k")))
        H = sum([psi(t)*f for psi, f in zip(self.psis, fs)])-(alpha1+alpha2*sp.Pow(u, 2))

        if self.phi:
            self.phi = sp.lambdify([x(t) for x in xs], self.phi)

        if H.diff(u, 2) == 0:
            if self.is_limit:
                if (alpha1 + alpha2 * sp.Pow(u, 2)) == 1:
                    u_opt = sp.sign(H.diff(u))
                else:
                    u_opt = sp.Piecewise((-1, H.diff(u)< -1), (1, H.diff(u)>1), (H.diff(u), True))
            else:
                u_opt = H.diff(u)
        else:
            if self.is_limit:
                u_opt_part = sp.solve(H.diff(u), u)[0]
                u_opt = sp.Piecewise((u_opt_part, sp.Abs(u_opt_part) <= 1), (1, u_opt_part > 1), (-1, u_opt_part < -1))
                # u_opt = sp.sign(H.diff(u))
            else:
                u_opt = sp.solve(H.diff(u), u)[0]

        if H.diff(u, 2) == 0:
            if self.is_limit:
                c = sp.Symbol('c')

                F =  sp.lambdify((*[x(t) for x in xs], *[psi(t) for psi in self.psis]),
                     [H.subs(u, c).diff(psi(t)).subs(c, u_opt) for psi in self.psis] + [-H.subs(u, c).diff(x(t)).subs(c, u_opt) for x in xs],'scipy',dummify=True)
            else:
                F = sp.lambdify((*[x(t) for x in xs], *[psi(t) for psi in self.psis]),
                                [H.subs(u, u_opt).diff(psi(t)) for psi in self.psis] + [-H.subs(u, u_opt).diff(x(t)) for
                                                                                        x in xs], ('scipy', 'numpy'),
                                dummify=True)
        else:
            F =  sp.lambdify((*[x(t) for x in xs], *[psi(t) for psi in self.psis]),
                 [H.subs(u, u_opt).diff(psi(t)) for psi in self.psis]+[-H.subs(u, u_opt).diff(x(t)) for x in xs],('scipy','numpy'),dummify=True)

        self.Hf = sp.lambdify((*[x(t) for x in xs], *[psi(t) for psi in self.psis]),H.subs(u, u_opt))
        def F_(*x):
            t, x0 = x
            return F(*x0)
    def solve(self):
        def Nev(psi0):
            # Проверка корректности конечного времени
            if psi0[-1] < 0:
                return [1e9] * (self.n+1)  # Возвращаем большое число, если время отрицательное

            # Интегрирование системы
            tn = np.linspace(t0, psi0[-1], self.num)
            sol = solve_ivp(F_, [t0, psi0[-1]], np.concatenate((b0,psi0[:-1])), t_eval=tn, method=self.diffmethod, atol=self.eps)
            psiT = sol.y[:, -1]  # Последнее значение решения
            # Выделение переменных состояния
            xf = psiT[:self.n]
            # Возвращаем разность между текущим и целевым состоянием
            if self.isBkFun:
                if self.tkCurrent:
                    return np.concatenate(([self.phi(*xf)]*self.n,[self.Hf(*xf,*psiT[self.n:-1], args[-1])]))
                else:
                    return np.concatenate(([self.phi(*xf)]*self.n, [self.Hf(*xf, *psiT[self.n:])]))
            else:
                if self.tkCurrent:
                    return np.concatenate(([i - k for i, k in zip(xf, self.bk)],[self.Hf(*self.bk,*psiT[self.n:-1], args[-1])]))
                else:
                    return np.concatenate(([i - k for i, k in zip(xf, self.bk)], [self.Hf(*self.bk, *psiT[self.n:])]))

        psi0 = args
        if self.tkCurrent:
            self.psi0 = args[:-1]
        else:
            self.psi0 = args
        y0 = root(Nev, psi0, method=self.method, tol=self.eps)['x']
        u = sp.lambdify((*[x(t) for x in xs], *[psi(t) for psi in self.psis]), u_opt, modules='numpy')
        if self.tkCurrent:
            tn = np.linspace(t0, args[-1], self.num)
        else:
            tn = np.linspace(t0, y0[-1], self.num)


        # if self.findBalancePoint:
        #     b = sp.solve(fs, [x(t) for x in xs])

        #     c = sp.symbols(f'c_1:{self.n+1}')
        #     self.bp = {sp.Symbol(x.name):ci for x,ci in zip(xs,c)}
        #     defaultbp = self.bp
        #     su = sp.Function('u')
        #     xsym = {x(t):sp.Symbol(x.name) for x in xs}

        #     self.bp.update(sp.solve([sp.Eq(f.subs(xsym).subs({su(t):sp.Symbol(su.name)}),0) for f in fs]))

        #     kwargs["balanceLine"] = False
        #     kwargs["balancePoint"] = False
        #     for key, value in self.bp.items():
        #         self.bp[key] = self.bp[key].subs({ci:0 for ci in c})
        #         if self.bp[key] in c:
        #             kwargs["balanceLine"] = True
        #     if not kwargs["balanceLine"]:
        #         kwargs["balancePoint"] = True
        #     self.bp = list(self.bp.values())[:-1]
        #     kwargs["balanceFunction"] = sp.lambdify(c, self.bp)


        y_toch = solve_ivp(F_, [t0, tn[-1]], np.concatenate((b0, y0[:-1])), t_eval=tn, method = self.diffmethod).y.T
        self.psin = y_toch[:,self.n: self.n*2]

        xn = y_toch[:,0:self.n]
        self.xk = xn[-1]
        
        if self.isBkFun:
            if self.phi(*self.xk)>=0.1:
                return
        else:
            if max([abs(i-j)**2 for i,j in zip(self.xk,self.bk)])>=0.1:
                return
            

        u = [u(*x, *psi) for x, psi in zip(xn, self.psin)]


        if self.tkCurrent:
            self.tk = args[-1]
        else:
            self.tk = y0[-1]