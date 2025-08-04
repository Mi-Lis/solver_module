import logging
from .tools import logged
class Solver:
    __res = {}
    def solve(self, *args): 
        self.prepareData()
        self.findMinimum()
        self.findOptimal()
        self.checkResult()
        self.updateResult()

    @logged
    def prepareData(self):
        logging.info(f"Start {self.prepareData.__name__}")
        ...
    @logged
    def findMinimum(self):
        logging.info(f"Start {self.findMinimum.__name__}")
        ...
    @logged
    def findOptimal(self):
        logging.info(f"Start {self.findOptimal.__name__}")
        ...
    @logged
    def checkResult(self):
        logging.info(f"Start {self.checkResult.__name__}")
        ...
    @logged
    def updateResult(self):
        logging.info(f"Start {self.updateResult.__name__}")
        ...
    @property
    def result(self):
        ...
    @property
    def name(self):
        return self.__name__
class Result:
    def get(solver:Solver):
        return solver.result
