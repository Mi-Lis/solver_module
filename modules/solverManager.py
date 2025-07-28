import logging


from solvers import _solvers

from solvers.solver import Solver



class SolverManager:
    typeAlg:str
    def __init__(self, typeAlg, fs, bs, psi0=[], **kwargs):
        self.typeAlg = typeAlg
        self.fs = fs
        self.bs = bs
        self.psi0 = psi0
        self.kwargs = kwargs
        pass

    def get_solver(self) -> Solver:
        for solver in _solvers:
            if solver.__qualname__ == self.typeAlg:
                return solver(self.fs, self.bs, self.psi0, **self.kwargs)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, filename="solver.log",filemode="w",
                    format="%(asctime)s %(levelname)s %(message)s")

    solver = SolverManager("BellmanCompute", [0,0], [0,0], [0,0]).get_solver()
    print(type(solver.solve([])))