import numpy as np
import sympy as sp
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from sympy.abc import t
from functools import reduce
from itertools import product
import operator


class BellmanCompute:
    def __init__(self, dfs, fs, bs,args=(), **kwargs):
        self.eps = 1e-5
        self.num = 50
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
        lam, tk = args
        self.psi0 = args
        t0 = bs[-1]
        xs = [sp.Function(f'x_{i}') for i in range(1, self.n+1)]
        bordersk = bs[0:self.n]
        S = sp.Function('S')(t, *[x(t) for x in xs])
        u  = sp.Symbol('u')
        fys = fs[:self.n]
        fys = [f.subs({x(t):sp.Symbol(x.name) for x in xs}).subs(sp.Function('u')(t), u) for f in fys]
        ds = [sp.Derivative(S, x(t)) for x in xs]
        self.functional = sp.Integral(u**2, (t, 0, sp.Symbol("t_k")))+sp.Symbol("lambda")*fs[-1]
        H = sum(f*s for f, s in zip(fys, ds))+u**2
        if self.findBalancePoint:
            b = sp.solve(fs, [x(t) for x in xs])
            c = sp.symbols(f'c_1:{self.n+1}')
            self.bp = {sp.Symbol(x.name):ci for x,ci in zip(xs,c)}
            defaultbp = self.bp
            su = sp.Function('u')
            xsym = {x(t):sp.Symbol(x.name) for x in xs}

            self.bp.update(sp.solve([sp.Eq(f.subs(xsym).subs({su(t):sp.Symbol(su.name)}),0) for f in fys]))
            kwargs["balanceLine"] = False
            kwargs["balancePoint"] = False
            for key, value in self.bp.items():
                self.bp[key] = self.bp[key].subs({ci:0 for ci in c})
                if self.bp[key] in c:
                    kwargs["balanceLine"] = True
            if not kwargs["balanceLine"]:
                kwargs["balancePoint"] = True

            self.bp = list(self.bp.values())[:-1]
            kwargs["balanceFunction"] = sp.lambdify(c, self.bp)
        power_sym = sp.symbols(f'p_1:{len(xs)+1}')
        xij = []
        c = 0
        ss = []
        for ij in product(*(range(len(fys)+1) for _ in range(len(xs)))):
            if sum(ij) <= len(fys):
                xij.append(reduce(operator.mul, [sp.Symbol(x.name)**powsym for x, powsym in zip(xs, power_sym)]))
                xij[-1] = xij[-1].subs({p:i for p, i in zip(power_sym, ij)})
                c+=1
                ss.append(sp.Function(f's_{c}'))

        expr = sum(s(t)*x for s, x in zip(ss, xij))

        u_kr = sp.solve(H.diff(u), u)[0]

        dxexpr = {S.diff(x(t)):expr.diff(sp.Symbol(x.name)) for x in xs}
        H = H.subs(u, u_kr).subs(dxexpr).expand()
        ST = lam*fs[-1].subs(t, tk).subs({x(tk):sp.Symbol(x.name) for x in xs})
        expandExpr = expr.subs({s(t):sp.Symbol(s.name) for s in ss})-ST
        systemEq = []
        for x in xij[::-1]:
            systemEq.append(expandExpr.coeff(x))
            expandExpr -= expandExpr.coeff(x) * x
            expandExpr = expandExpr.expand()

        y0Solve = sp.solve(systemEq, [sp.Symbol(s.name) for s in ss],dict=True)[0]
        y0 = []
        for key, value in y0Solve.items():
            y0.append(y0Solve[key])

        PolyH = sp.Poly(H, [sp.Symbol(x.name) for x in xs])

        ns = [sp.lambdify((*(sp.Symbol(s.name) for s in ss),t), f.subs({s(t):sp.Symbol(s.name) for s in ss})) for f in list(PolyH.as_dict().values())]
        def _ns(*x):
            t,y = x
            return [n(*y,t) for n in ns]

        nt = np.linspace(t0, tk, self.num)

        f = solve_ivp(_ns,
                      [t0, tk],
                      y0,
                      t_eval=nt,
                      method=self.diffmethod,
                      rtol=1e-8,  # Уменьшаем допустимую ошибку
                      atol=1e-8,
                      ).y.T

        sSpline = []
        for i in range(len(f[0])):
            sSpline.append(CubicSpline(nt, f[:,i]))
        for i, f in enumerate(fys):
            fys[i] = f.subs(u, u_kr.subs(dxexpr)).subs({sp.Symbol(x.name):x(t) for x in xs})
        lfys = []

        for f in fys:
            lfys.append(sp.lambdify((t, *(s(t) for s in ss), *(x(t) for x in xs)), f.subs(dxexpr)))
        def dfys(t, x):
            return [f(t, *(s(t) for s in sSpline), *x) for f in lfys]
        u = sp.lambdify((t, *(s(t) for s in ss), *(sp.Symbol(x.name) for x in xs)), u_kr.subs(dxexpr))
        xs = solve_ivp(dfys, [t0, tk], bordersk, t_eval=nt, method=self.diffmethod).y.T
        _u = u(nt, *[s(nt) for s in sSpline], *[xs[0:,i] for i in range(self.n)])

        bbs = ([0]*self.n+bs[0:self.n])
        bbs.reverse()

        self.u_opt = sp.latex(sp.Eq(sp.Function('u')(sp.Symbol('t')),u_kr.subs(dxexpr)))
        self.H = sp.latex(H)
