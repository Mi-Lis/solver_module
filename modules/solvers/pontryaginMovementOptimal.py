
from pontryaginBase import PontryaginCompute


class PontryaginMovementOptimal(PontryaginCompute):
    def __init__(self, fs, bs, args = (), **kwargs):
        alpha1 = 1
        alpha2 = 0
        super().__init__(fs, bs, alpha1, alpha2, args, **kwargs)