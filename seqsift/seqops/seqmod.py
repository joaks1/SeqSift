#! /usr/bin/env python

import sys
import os
from string import maketrans

from seqsift.seqops.sequtils import copy_seq_metadata
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def remove_gaps(seq_iter, gap_characters=['-']):
    gap_chars = ''.join(gap_characters)
    return seq_mod(seq_iter, del_chars=gap_chars)

def seq_mod(seq_iter,
        from_chars='',
        to_chars='',
        del_chars=''):
    '''
    Modify sequences. Each sequence in `seq_iter` will have the characters
    in the `from_chars` string mapped to the characters in the `to_chars`
    string, and any characters in the `del_chars` string will be removed.
    '''
    if len(from_chars) != len(to_chars):
        raise ValueError('from and to characters must have same length')
    table = None
    if len(from_chars) > 0:
        table = maketrans(from_chars, to_chars)
    for seq in seq_iter:
        yield copy_seq_metadata(seq,
                new_seq=str(seq.seq).translate(table, del_chars))

