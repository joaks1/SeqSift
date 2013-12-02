#! /usr/bin/env python

import sys
import os

from seqsift.align import global_align
from seqsift.utils.dataio import BufferedIter
from seqsift.utils.errors import AlignmentError
from seqsift.utils.alphabets import DnaAlphabet
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def column_frequencies(seq_iter, character_list=['-','?']):
    seqs = BufferedIter(seq_iter)
    char_list = [c.lower() for c in character_list]
    column_counts = []
    align_length = None
    for i, seq_record in enumerate(seqs):
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
    
def distance(seq1, seq2,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        similarity_matrix = {
            ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
        },
        gap_cost = -5):
    diffs, l  = get_differences(seq1 = seq1, seq2 = seq2,
            aligned = aligned,
            ignore_gaps = ignore_gaps,
            alphabet = alphabet,
            similarity_matrix = similarity_matrix,
            gap_cost = gap_cost)
    return len(diffs)

def per_site_distance(seq1, seq2,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        similarity_matrix = {
            ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
        },
        gap_cost = -5):
    diffs, l  = get_differences(seq1 = seq1, seq2 = seq2,
            aligned = aligned,
            ignore_gaps = ignore_gaps,
            alphabet = alphabet,
            similarity_matrix = similarity_matrix,
            gap_cost = gap_cost)
    return len(diffs) / float(l)

def get_differences(seq1, seq2,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        similarity_matrix = {
            ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
        },
        gap_cost = -5):
    if not alphabet:
        alphabet = DnaAlphabet()
    if not aligned:
        seq1, seq2 = global_align(
                [x for x in seq1 if x != alphabet.gap],
                [x for x in seq2 if x != alphabet.gap],
                similarity_matrix = similarity_matrix,
                gap_cost = gap_cost)
    if len(seq1) != len(seq2):
        raise AlignmentError('Sequences are not aligned')
    residue_codes = alphabet.all_residue_codes
    diffs = {}
    for i in range(len(seq1)):
        if seq1[i].upper() == seq2[i].upper():
            continue
        if (seq1[i] == alphabet.gap) or (seq2[i] == alphabet.gap):
            if not ignore_gaps:
                diffs[i] = (seq1[i].upper(), seq2[i].upper())
            continue
        s1 = set(residue_codes[seq1[i].upper()])
        s2 = set(residue_codes[seq2[i].upper()])
        if not s1.intersection(s2):
            diffs[i] = (seq1[i].upper(), seq2[i].upper())
    return diffs, len(seq1)

