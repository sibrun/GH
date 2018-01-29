import cProfile
import functools
import os
import StoreLoad as SL

reload(SL)


def profile(dir_name):
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            path = os.path.join(dir_name, func.__name__ + ".profile")
            SL.generate_path(path)
            prof = cProfile.Profile()
            re = prof.runcall(func, *args, **kwargs)
            prof.dump_stats(path)
            return re
        return wrapper
    return inner
