import logging

from models.formula import Task


from .solvers import _solvers

from .solvers.solver import Solver



class SolverManager:
    typeAlg:str
    def __init__(self):

        pass

    def get_solver(self, task:Task, **kwargs) -> Solver:
        # self.typeAlg = typeAlg
        # self.fs = fs
        # self.bs = bs
        # self.psi0 = psi0
        self.kwargs = kwargs
        for solver in _solvers:
            if solver.__qualname__ == task.type:
                return solver([task.equation1,
                               task.equation2], 
                               task.border_boundary, 
                               task.psi_boundary+[task.t_end],
                               **self.kwargs)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, filename="solver.log",filemode="w",
                    format="%(asctime)s %(levelname)s %(message)s")

    solver = SolverManager("BellmanCompute", [0,0], [0,0], [0,0]).get_solver()
    print(type(solver.solve([])))