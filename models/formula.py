from typing import List

from pydantic import BaseModel

# class Formula(BaseModel):
#     latex:str

class Task(BaseModel):
    type:str="PontryaginCompute"
    equation1:str
    equation2:str
    border_boundary:List[float]
    psi_boundary:List[float]
    t_end:float
