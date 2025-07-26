from pontryaginBase import PontryaginCompute

class PontryaginCustomOptimal(PontryaginCompute):
    def __init__(self, fs, bs:list, args=(), **kwargs):
        alpha2 = bs[-2]
        alpha1 = bs[-1]
        bs = bs[:-2]
        print(bs)
        super().__init__(fs, bs, alpha1, alpha2, args, **kwargs)