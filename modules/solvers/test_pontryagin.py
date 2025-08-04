from pontryaginBase import PontryaginCompute
import sympy as sp



def test_simple():
    t = sp.Symbol('t')
    x2 = sp.Function('x_2')
    u = sp.Function('u')(t)
    fs = [x2(t), x2(t)-u, u**2]
    bs = [1, 0, 0, 0, 0]
    args = [3,5, 3]
    pontryagin = PontryaginCompute(fs, bs, args)
    pontryagin.solve()
    print(pontryagin.result)




if __name__ == '__main__':
    test_simple()