from pydantic import BaseModel

class Task(BaseModel):
    id:int
    typeTask:str
    functions:list[dict]
    boundary:list[list[float]]
    psi0:list[float]