#! /usr/bin/env python

import sys
import os
import logging

_LOGGING_LEVEL_ENV_VAR = "SEQSIFT_LOGGING_LEVEL"

def get_logging_level():
    if _LOGGING_LEVEL_ENV_VAR in os.environ:
        if os.environ[_LOGGING_LEVEL_ENV_VAR].upper() == "NOTSET":
            logging_level = logging.NOTSET
        elif os.environ[_LOGGING_LEVEL_ENV_VAR].upper() == "DEBUG":
            logging_level = logging.DEBUG
        elif os.environ[_LOGGING_LEVEL_ENV_VAR].upper() == "INFO":
            logging_level = logging.INFO
        elif os.environ[_LOGGING_LEVEL_ENV_VAR].upper() == "WARNING":
            logging_level = logging.WARNING
        elif os.environ[_LOGGING_LEVEL_ENV_VAR].upper() == "ERROR":
            logging_level = logging.ERROR
        elif os.environ[_LOGGING_LEVEL_ENV_VAR].upper() == "CRITICAL":
            logging_level = logging.CRITICAL
        else:
            logging_level = logging.WARNING
    else:
        logging_level = logging.WARNING
    return logging_level

def get_logger(name = 'seqsift'):
    log = logging.getLogger(name)
    level = get_logging_level()
    log.setLevel(level)
    h = logging.StreamHandler()
    h.setLevel(level)
    log.addHandler(h)
    return log
