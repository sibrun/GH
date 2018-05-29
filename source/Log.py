"""Provide a common logger."""

import logging
import os
import StoreLoad as SL
import Parameters


logger = logging.getLogger('graph_homology')
logger.addHandler(logging.NullHandler())


def set_log_file(log_file):
    log_path = os.path.join(Parameters.log_dir, log_file)
    SL.generate_path(log_path)
    logger.addHandler(logging.FileHandler(log_path))


log_levels_dict = {'info': logging.INFO, 'warning': logging.WARNING, 'error': logging.ERROR}


def set_log_level(level):
    log_level = log_levels_dict.get(level)
    if log_level is None:
        raise ValueError('not supported log level')
    else:
        logger.setLevel(log_level)