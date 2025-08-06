from itertools import product
from typing import List
import logging

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root
from sympy.parsing.latex import parse_latex_lark
from lark import Tree

from .solver import Solver
from .tools import logged


class PontryaginCompute(Solver):
    __name__ = "Метод Понтрягина"
    __res = {
        "tk":None,
        "xn":None,
        "un":None
    }
    def __init__(self, fs:List[sp.Equality], bs:List[float], args=(),u_opt=None, **kwargs) -> None:
        super().__init__()
        self.eps = 1e-5
        u = sp.Symbol('u')
        self.t = sp.Symbol('t')
        self.num = 50
        self.n = 2
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
        self.fs = [parse_latex_lark(f) for f in fs]
        print(type(self.fs))
        self.functional = 1
    
    def get_value_from_tree_lark(self, tree):
        value = []
        for child in tree.children:
                if isinstance(child, Tree) and child.data == '_ambig':
                   value.append(child.children[0])
                else:
                    value.append(child)
        return value[0]
    
    @logged
    def prepareData(self):
        self.xs = [sp.Symbol(r"x"), sp.Symbol(r"y")]
        self.psis = list(sp.symbols(f'psi1:{self.n+1}', cls=sp.Function))
        t = sp.Symbol('t')
        u = sp.Symbol('u')
        H = sum([psi(t)*f for psi, f in zip(self.psis, self.fs)])-(self.functional)
        print(H)

        if self.phi:
            self.phi = sp.lambdify([x(t) for x in self.xs], self.phi)

        if H.diff(u, 2) == 0:
                # u_opt = sp.Piecewise((-1, H.diff(u)< -1), (1, H.diff(u)>1), (H.diff(u), True))
                self.u_opt = H.diff(u)
        else:
                self.u_opt = sp.solve(H.diff(u), u)[0]
        self.F =  sp.lambdify((*[_x for _x in self.xs], *[psi(t) for psi in self.psis]),
                 [H.subs(u, self.u_opt).diff(psi(t)) for psi in self.psis]+[-H.subs(u, self.u_opt).diff(x) for x in self.xs],dummify=True)

        self.Hf = sp.lambdify((*[_x for _x in self.xs], *[psi(t) for psi in self.psis]),H.subs(u, self.u_opt))
    
    def F_(self, *x):
            t, x0 = x
            print(x0)
            return self.F(*x0)
    @logged
    def findMinimum(self):

        def Nev(psi0):
            # Проверка корректности конечного времени
            if psi0[-1] < 0:
                return [1e9] * (self.n+1)  # Возвращаем большое число, если время отрицательное

            # Интегрирование системы
            tn = np.linspace(self.t0, psi0[-1], self.num)
            sol = solve_ivp(fun=self.F_, 
                            t_span=[self.t0, psi0[-1]], 
                            y0=np.concatenate((self.b0,psi0[:-1])), 
                            t_eval=tn, 
                            method=self.diffmethod, 
                            atol=self.eps)
            psiT = sol.y[:, -1]  # Последнее значение решения
            # Выделение переменных состояния
            xf = psiT[:self.n]
            # Возвращаем разность между текущим и целевым состоянием
            return np.concatenate(([i - k for i, k in zip(xf, self.bk)], [self.Hf(*self.bk, *psiT[self.n:])]))

        psi0 = self.psi0
        
        self.y0 = root(Nev, psi0, method=self.method, tol=self.eps)['x']
    
    @logged
    def findOptimal(self):
        self.tn = np.linspace(self.t0, self.y0[-1], self.num)
        y_toch = solve_ivp(self.F_, [self.t0, self.tn[-1]], np.concatenate((self.b0, self.y0[:-1])), t_eval=self.tn, method = self.diffmethod).y.T
        self.psin = y_toch[:,self.n:self.n*2]
        print(self.psin, self.psis)
        self.xn = y_toch[:,0:self.n]
        self.xk = self.xn[-1]

    @logged
    def checkResult(self):
            assert max([abs(i-j)**2 for i,j in zip(self.xk,self.bk)])>=0.1

    @logged 
    def updateResult(self):
        self.u = sp.lambdify((*[x for x in self.xs], *[psi(self.t) for psi in self.psis]), self.u_opt, modules='numpy')
        u = np.array([self.u(*x, *psi) for x, psi in zip(self.xn, self.psin)])


        self.__res["tk"] = self.tn
        self.__res["xn"] = self.xn
        self.__res["un"] = u

    @property
    def result(self):
        return self.__res