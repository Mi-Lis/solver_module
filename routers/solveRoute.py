from fastapi import APIRouter

from models.task import Task

from modules.solverManager import SolverManager

__all__ = ['solveRouter']

solveRouter = APIRouter("solve/")
solverManager = SolverManager()

@solveRouter.post("")
async def postTask(taskData:Task):
    await solverManager.create_task(taskData)
    ...


@solveRouter.get("{taskId}/")
async def getAnswer(taskId):
    await solverManager.solve(taskId)
    ...