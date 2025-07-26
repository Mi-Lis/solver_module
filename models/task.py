from pydantic import BaseModel

class Task(BaseModel):
    typeTask:str
    functions:list[dict]
    boundary:list[list[float]]