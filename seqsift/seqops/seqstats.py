#! /usr/bin/env python

import sys
import os

from seqsift.utils.dataio import BufferedIter
from seqsift.utils.errors import AlignmentError
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def column_frequencies(seq_iter, character_list=['-','?']):
    seqs = BufferedIter(seq_iter)
    char_list = [c.lower() for c in character_list]
    column_counts = []
    align_length = None
    for i, seq_record in enumerate(seqs.iter()):
        if align_length == None:
            align_length = len(seq_record)
            column_counts = [0] * align_length
        else:
            if len(seq_record) != align_length:
                raise AlignmentError('Sequence {0} has unexpected '
                        'length {1}.'.format(seq_record.name, len(seq_record)))
        for j, character in enumerate(seq_record.seq):
            if character.lower() in char_list:
                column_counts[j] += 1
    column_freqs = [count / float(i + 1) for count in column_counts]
    return column_freqs, seqs
    

