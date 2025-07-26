
from enum import Enum
from solvers import BellmanCompute, PontryaginCustomOptimal, PontryaginEnergyOptimal, PontryaginMovementOptimal


class SolverManager:
    class TASKS(Enum):
        naiskDvi = u"Наискорейшее движение точки в сопротивляющейся среде"
        linSys = u"Линейные системы с квадратичным критерием качества"
        enOpt = u"Энергетически оптимальное движение точки в сопротивляющейся среде"
        metodDyn = u"Метод динамического программирования"
    pass
    typeAlg:str
    def __init__(self, typeAlg, fs, bs, psi0=[], **kwargs):
        self.typeAlg = typeAlg
        self.fs = fs
        self.bs = bs
        self.psi0 = psi0
        self.kwargs = kwargs
        pass

    def solve(self):
        match self.typeAlg:
            case self.TASKS.naiskDvi.value:
                return PontryaginMovementOptimal(self.fs, self.bs, self.psi0, **self.kwargs).solve()
            case self.TASKS.linSys.value:
                return PontryaginCustomOptimal(self.fs, self.bs, self.psi0, **self.kwargs).solve()
            case self.TASKS.enOpt.value:
                return PontryaginEnergyOptimal(self.fs, self.bs, self.psi0, **self.kwargs).solve()
            case self.TASKS.metodDyn.value:
                return  BellmanCompute(self.fs, self.bs, self.psi0, **self.kwargs).solve()

