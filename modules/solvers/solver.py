import logging
from .tools import logged
class Solver:
    __res = {}
    def solve(self, *args): 

        try:
            self.prepareData()
        except:
            logging.error("Failed prepare data", exc_info=True)
        
        try:
            self.findMinimum()
        except: 
            logging.error("Failed find minimum", exc_info=True)
        

        try:
            self.findOptimal()
        except:
            logging.error("Failed find optimal", exc_info=True)
        
        try: 
            self.checkResult()
        except: 
            logging.error("Failed result", exc_info=True)
        
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
