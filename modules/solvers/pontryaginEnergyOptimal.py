from pontryaginBase import PontryaginCompute


class PontryaginEnergyOptimal(PontryaginCompute):
    def __init__(self, fs, bs, args = (), **kwargs):
        alpha1 = 1
        alpha2 = 1
        super().__init__(fs, bs, alpha1, alpha2, args, **kwargs)