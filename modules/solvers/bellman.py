import numpy as np
import sympy as sp
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from sympy.abc import t
from functools import reduce
from itertools import product
import operator
import logging
from .solver import Solver
from .tools import logged

class BellmanCompute(Solver):
    __name__ = "Метод Беллмана"
    def __init__(self, fs=None, bs=[0,0],args=(0,0), **kwargs):
        logging.info("init")
        self.u  = sp.Symbol('u')
        self.fs = fs
        self.eps = 1e-5
        self.num = 50
        self.n = 2
        self.method = 'hybr'
        self.diffmethod='RK45'
        self.findBalancePoint = False
        for key, value in kwargs.items():
            match key:
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
                    if value == 'True':
                        value = True
                    else:
                        value = False
                    if bool(value):
                        self.findBalancePoint=True
        self.lam, self.tk = args
        self.psi0 = args
        self.bordersk = bs[0:self.n]
        self.t0 = bs[-1]
        self.functionalCore = fs[-1]
        self.xs = [sp.Function(f'x_{i}') for i in range(1, self.n+1)]
   
    
    @logged
    def solve(self, *args):
        super().solve(*args)
    
    @logged
    def prepareData(self):
        super().prepareData()
        self.fys = [f.subs({x(t):sp.Symbol(x.name) for x in self.xs}).subs(sp.Function('u')(t), self.u) for f in self.fs[:self.n]]
        S = sp.Function('S')(t, *[x(t) for x in self.xs])

        
        ds = [sp.Derivative(S, x(t)) for x in self.xs]
        self.functional = sp.Integral(self.u**2, (t, 0, sp.Symbol("t_k")))+sp.Symbol("lambda")*self.functionalCore
        self.H = sum(f*s for f, s in zip(self.fys, ds))+self.u**2
        power_sym = sp.symbols(f'p_1:{self.n+1}')
        xij = []
        c = 0
        ss = []
        for ij in product(*(range(self.n+1) for _ in range(self.n))):
            if sum(ij) <= self.n:
                xij.append(reduce(operator.mul, [sp.Symbol(x.name)**powsym for x, powsym in zip(self.xs, power_sym)]))
                xij[-1] = xij[-1].subs({p:i for p, i in zip(power_sym, ij)})
                c+=1
                ss.append(sp.Function(f's_{c}'))

        expr = sum(s(t)*x for s, x in zip(ss, xij))

        self.u_kr = sp.solve(self.H.diff(self.u), self.u)[0]

        self.dxexpr = {S.diff(x(t)):expr.diff(sp.Symbol(x.name)) for x in self.xs}
        self.H = self.H.subs(self.u, self.u_kr).subs(self.dxexpr).expand()
        ST = self.lam*self.functionalCore.subs(t, self.tk).subs({x(self.tk):sp.Symbol(x.name) for x in self.xs})
        expandExpr = expr.subs({s(t):sp.Symbol(s.name) for s in ss})-ST
        systemEq = []
        for x in xij[::-1]:
            systemEq.append(expandExpr.coeff(x))
            expandExpr -= expandExpr.coeff(x) * x
            expandExpr = expandExpr.expand()

        y0Solve = sp.solve(systemEq, [sp.Symbol(s.name) for s in ss],dict=True)[0]
        self.y0 = []
        for key, value in y0Solve.items():
            self.y0.append(y0Solve[key])

        PolyH = sp.Poly(self.H, [sp.Symbol(x.name) for x in self.xs])

        self.ns = [sp.lambdify((*(sp.Symbol(s.name) for s in ss),t), f.subs({s(t):sp.Symbol(s.name) for s in ss})) for f in list(PolyH.as_dict().values())]

    @logged
    def findMinimum(self):
        def _ns(*x):
            t,y = x
            return [n(*y,t) for n in self.ns]
        self.nt = np.linspace(self.t0, self.tk, self.num)

        self.f = solve_ivp(_ns,
                      [self.t0, self.tk],
                      self.y0,
                      t_eval=self.nt,
                      method=self.diffmethod,
                      rtol=1e-8,  # Уменьшаем допустимую ошибку
                      atol=1e-8,
                      ).y.T


    @logged
    async def findOptimal(self):
        sSpline = []
        for i in range(len(f[0])):
            sSpline.append(CubicSpline(self.nt, self.f[:,i]))
        for i, f in enumerate(self.fys):
            self.fys[i] = self.f.subs(u, self.u_kr.subs(self.dxexpr)).subs({sp.Symbol(x.name):x(t) for x in self.xs})
        lfys = []

        for f in self.fys:
            lfys.append(sp.lambdify((t, *(s(t) for s in self.ss), *(x(t) for x in self.xs)), f.subs(self.dxexpr)))

        u = sp.lambdify((t, *(s(t) for s in self.ss), *(sp.Symbol(x.name) for x in self.xs)), self.u_kr.subs(self.dxexpr))
        def dfys(t, x):
            return [f(t, *(s(t) for s in sSpline), *x) for f in lfys]
        self.xs = await solve_ivp(dfys, [self.t0, self.tk], self.bordersk, t_eval=self.nt, method=self.diffmethod).y.T
    @logged
    def checkResult(self):
        return
    @logged
    def updateResult(self):
        self._u = self.u(self.nt, *[s(self.nt) for s in self.sSpline], *[self.xs[0:,i] for i in range(self.n)])

        bbs = ([0]*self.n+self.bs[0:self.n])
        bbs.reverse()
        self.u_opt = sp.latex(sp.Eq(sp.Function('u')(sp.Symbol('t')),self.u_kr.subs(self.dxexpr)))
        self.H = sp.latex(self.H)
        self.__res["H"] = self.H
        self.__res["u_opt"] = self.u_opt
    @property
    def result(self):
        return self.__res