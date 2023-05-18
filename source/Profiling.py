"""Provides a conditional decorator to profile functions."""

import cProfile
import functools
import os
import StoreLoad as SL


def cond_decorator(flag, dec):
    def decorate(func):
        return dec(func) if flag else func
    return decorate


def profile(dir_name):
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            prof = cProfile.Profile()
            re = prof.runcall(func, *args, **kwargs)
            path = os.path.join(dir_name, func.__name__ + ".prof")
            SL.generate_path(path)
            prof.dump_stats(path)
            return re
        return wrapper
    return inner
