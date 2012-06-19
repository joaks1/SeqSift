#! /usr/bin/env python

import sys
import os

from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def length_filter(seq_iter, min_length=0, max_length=None):
    if max_length and (max_length < min_length):
        raise ValueError("max_length cannot be smaller than min_length")
    for seq in seq_iter:
        if max_length:
            if len(seq) <= max_length and len(seq) >= min_length:
                yield seq
        else:
            if len(seq) >= min_length:
                yield seq
