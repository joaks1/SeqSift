#! /usr/bin/env python

import sys
import os
import random
import errno
import subprocess
import string

from seqsift.utils import GLOBAL_RNG

def random_str(length=8,
        char_pool=string.ascii_letters + string.digits):
    return ''.join(random.choice(char_pool) for i in range(length))

def sample_iter(iterable, sample_size, exclude = [], exclude_attribute = None,
        rng = None):

    def get_attribute(el, attribute=None):
        if attribute is None:
            return el
        return getattr(el, attribute)

    exclude = [get_attribute(x, exclude_attribute) for x in exclude]
    if not rng:
        rng = GLOBAL_RNG
    samples = []
    iterator = iter(iterable)
    try:
        for i in range(sample_size):
            s = next(iterator)
            while get_attribute(s, exclude_attribute) in exclude:
                s = next(iterator)
            samples.append(s)
    except StopIteration:
        raise ValueError("Sample size {0} is larger than population".format(
                sample_size))
    rng.shuffle(samples)
    for idx, el in enumerate(iterator, sample_size):
        if get_attribute(el, exclude_attribute) in exclude:
            continue
        r = rng.randint(0, idx)
        if r < sample_size:
            samples[r] = el
    return samples

def mkdr(path):
    """
    Creates directory `path`, but suppresses error if `path` already exists.
    """
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise

def get_external_tool(exe_file):
    """
    Uses `subprocess.Popen` to check system for `exe_file`. If found,
    `exe_file` is returned, else `None` is returned.

    """
    try:
        p = subprocess.Popen([exe_file],
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE)
        p.terminate()
    except subprocess.CalledProcessError:
        return exe_file
    except OSError:
        raise Exception('Unable to find executable '
            '{0!r}'.format(exe_file))
    return exe_file

def get_tool_full_path(name):
    p = get_external_tool(name)
    if os.path.isfile(p):
        return p
    return which(name)

def get_new_path(path, max_attempts = 1000):
    path = os.path.abspath(os.path.expandvars(os.path.expanduser(path)))
    if not os.path.exists(path):
        f = open(path, 'w')
        f.close()
        return path
    attempt = 0
    while True:
        p = '-'.join([path, str(attempt)])
        if not os.path.exists(p):
            f = open(p, 'w')
            f.close()
            return p
        if attempt >= max_attempts:
            raise Exception('failed to get unique path')
        attempt += 1

def is_executable(path):
    return os.path.isfile(path) and bool(get_external_tool(path))

def which(exe):
    pth, name = os.path.split(exe)
    if pth:
        if is_executable(exe):
            return exe
    else:
        for p in os.environ['PATH'].split(os.pathsep):
            p = p.strip('"')
            exe_path = os.path.join(p, exe)
            if is_executable(exe_path):
                return exe_path
    return None

