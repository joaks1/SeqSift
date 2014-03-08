#! /usr/bin/env python

import os
import argparse

from seqsift.utils.fileio import expand_path
from seqsift.utils.functions import is_executable

def arg_is_path(path):
    try:
        if not os.path.exists(path):
            raise
    except:
        msg = 'path {0!r} does not exist'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_file(path):
    try:
        if not os.path.isfile(path):
            raise
    except:
        msg = '{0!r} is not a file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_dir(path):
    try:
        if not os.path.isdir(path):
            raise
    except:
        msg = '{0!r} is not a directory'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_executable(path):
    try:
        if not is_executable(path):
            raise
    except:
        msg = '{0!r} is not an executable'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_nonnegative_int(i):
    try:
        if int(i) < 0:
            raise
    except:
        msg = '{0!r} is not a non-negative integer'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return int(i)

def arg_is_positive_float(i):
    try:
        if float(i) <= 0.0:
            raise
    except:
        msg = '{0!r} is not a positive real number'.format(i)
        raise argparse.ArgumentTypeError(msg)
    return float(i)

