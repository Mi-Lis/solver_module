from pydantic import BaseModel

class Image(BaseModel):
    url: str
    name: str


class Item(BaseModel):
    name: str
    image: Image | None = None
