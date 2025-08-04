import logging
from functools import wraps


def logged(func):
    @wraps(func)
    def wrap(*args, **kwargs):
        logging.info(f"Started {func.__name__}")
        func(*args, **kwargs)
        # try: func(*args, **kwargs)
        # except:
        #     logging.error(f"Failed exec {func.__name__}")
    return wrap