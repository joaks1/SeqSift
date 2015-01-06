#! /usr/bin/env python

import sys
import os
import itertools

from seqsift.seqops import sequtils, seqstats
from seqsift.utils import dataio
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def duplicate_id_filter(seq_iter):
    ids = set()
    checked_seqs = dataio.BufferedIter()
    for seq in seq_iter:
        if seq.id in ids:
            for prev_seq in checked_seqs:
                if prev_seq.id == seq.id:
                    break
            if not sequtils.sequences_are_equal(prev_seq, seq):
                    raise Exception('Found {0} more than once, but with '
                            'different data!'.format(seq.id))
            continue
        else:
            checked_seqs.append(seq)
        ids.add(seq.id)
        yield seq

def id_filter(seq_iter, ids_to_exclude):
    for seq in seq_iter:
        if seq.id in ids_to_exclude:
            continue
        yield seq

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

def column_filter(seq_iter, character_list=['?','-'], max_frequency=1.0):
    col_freqs, seqs = seqstats.column_frequencies(seq_iter,
            character_list=character_list)
    cols_to_keep = [p < max_frequency for p in col_freqs]
    for seq in seqs:
        new_seq = itertools.compress(str(seq.seq), cols_to_keep)
        yield sequtils.copy_seq_metadata(seq, new_seq=''.join(new_seq))

def row_filter(seq_iter, character_list=['?','-'], max_frequency=1.0):
    char_list = [c.lower() for c in character_list]
    for seq in seq_iter:
        char_count = 0
        for char in str(seq.seq):
            if char.lower() in char_list:
                char_count += 1
        freq = char_count/float(len(seq))
        if freq < max_frequency:
            yield seq

