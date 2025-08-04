from typing import Annotated
from fastapi import APIRouter
from fastapi import Depends
from models.task import Task
from sqlalchemy.orm import Session

from modules.solverManager import SolverManager
from modules.taskManager import TaskManager


__all__ = ['solveRouter']

solveRouter = APIRouter(prefix="/solve/")

solverManager = SolverManager()
taskManager = TaskManager()

@solveRouter.post("")
async def postTask(taskData:Task):
    await taskManager.create_task(taskData)
    ...

def getTask(taskId):
    ...

@solveRouter.get("{taskId}/")
async def getAnswer(taskId, db:Annotated[Session, Depends(getTask)]):
    
    task = db.query(Task).filter(Task.id == taskId).first()
    solver = solverManager.get_solver(task.type, task.fs, task.bs)

    await solver.solve()
    
    db.commit()  
    db.refresh(task)
    
    return {"result":solver.result}