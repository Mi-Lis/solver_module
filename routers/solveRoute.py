from typing import Annotated
from fastapi import APIRouter, Request
from fastapi import Depends
from fastapi.responses import HTMLResponse
from models.formula import Task
# from sqlalchemy.orm import Session

from modules.ploters.plot import Plot
from modules.solverManager import SolverManager
from modules.taskManager import TaskManager


__all__ = ['solveRouter']

solveRouter = APIRouter(prefix="/solve")

solverManager = SolverManager()
taskManager = TaskManager()

@solveRouter.post("/", response_class=HTMLResponse)
async def postTask(task:Task):
    # return request
    # await taskManager.create_task(taskData)
    # return {"task":task}
    solver = solverManager.get_solver(task)
    solver.solve()
    plot = Plot()
    plot.post(**solver.result).savefig("plt.png")
    html_content = """<html>
        <head>
            <title></title>
        </head>
        <body>
        <img src="plt.png">
        <h1>Hello World</h1>
        </body>
    </html>
    """
    return HTMLResponse(content=html_content, status_code=200)


def getTask(taskId):
    ...

# @solveRouter.get("{taskId}/")
# async def getAnswer(taskId, db:Annotated[Session, Depends(getTask)]):
    
#     task = db.query(Task).filter(Task.id == taskId).first()
    

#     await solver.solve()
    
#     db.commit()  
#     db.refresh(task)
    
#     return {"result":solver.result}