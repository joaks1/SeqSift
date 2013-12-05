#! /usr/bin/env python

import sys
import os
import random

from seqsift.utils import GLOBAL_RNG

def sample_iter(iterable, sample_size, rng = None):
    if not rng:
        rng = GLOBAL_RNG
    samples = []
    iterator = iter(iterable)
    try:
        for i in xrange(sample_size):
            samples.append(iterator.next())
    except StopIteration:
        raise ValueError("Sample size {0} is larger than population".format(
                sample_size))
    rng.shuffle(samples)
    for idx, el in enumerate(iterator, sample_size):
        r = rng.randint(0, idx)
        if r < sample_size:
            samples[r] = el
    return samples
    
