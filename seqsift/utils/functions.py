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

