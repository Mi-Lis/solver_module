import logging
from functools import wraps


def logged(func):
    @wraps(func)
    def wrap(*args, **kwargs):
        logging.info(f"Started {func.__name__}")
        return func
    return wrap