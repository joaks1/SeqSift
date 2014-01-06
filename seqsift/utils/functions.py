#! /usr/bin/env python

import sys
import os
import random

from seqsift.utils import GLOBAL_RNG

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
        for i in xrange(sample_size):
            s = iterator.next()
            while get_attribute(s, exclude_attribute) in exclude:
                s = iterator.next()
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
    except OSError, e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise

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
    return os.path.isfile(path) and os.access(path, os.X_OK)

def which(exe):
    if is_executable(exe):
        return exe
    name = os.path.basename(exe)
    for p in os.environ['PATH'].split(os.pathsep):
        p = p.strip('"')
        exe_path = os.path.join(p, name)
        if is_executable(exe_path):
            return exe_path
    return None

