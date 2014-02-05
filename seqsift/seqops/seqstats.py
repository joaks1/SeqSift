#! /usr/bin/env python

import sys
import os
import itertools

from seqsift.align import align, align_pair
from seqsift.utils import functions
from seqsift.utils.dataio import BufferedIter
from seqsift.utils.errors import AlignmentError
from seqsift.utils.alphabets import DnaAlphabet
from seqsift.seqops import sequtils
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

def pairwise_distance_iter(seq_iter,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle']):
    seqs = BufferedIter(seq_iter)
    for seq1, seq2 in itertools.combinations(seqs, 2):
        d, drc = get_distances(
                seq1 = seq1,
                seq2 = seq2,
                per_site = per_site,
                aligned = aligned,
                ignore_gaps = ignore_gaps,
                alphabet = alphabet,
                aligner_tools = aligner_tools)
        yield seq1, seq2, d, drc

def sample_distance_iter(seq_iter,
        sample_size,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle'],
        rng = None):
    seqs = BufferedIter(seq_iter)
    seqs_to_sample = BufferedIter(seqs)
    for seq1 in seqs:
        samples = functions.sample_iter(iterable = seqs_to_sample,
                sample_size = sample_size,
                exclude = [seq1],
                exclude_attribute = 'id',
                rng = rng)
        for seq2 in samples:
            assert seq1.id != seq2.id
            d, drc = get_distances(
                    seq1 = seq1,
                    seq2 = seq2,
                    per_site = per_site,
                    aligned = aligned,
                    ignore_gaps = ignore_gaps,
                    alphabet = alphabet,
                    aligner_tools = aligner_tools)
            yield seq1, seq2, d, drc

def get_distances(seq1, seq2,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle']):
        d = distance(
                seq1 = seq1,
                seq2 = seq2,
                per_site = per_site,
                aligned = aligned,
                ignore_gaps = ignore_gaps,
                alphabet = alphabet,
                aligner_tools = aligner_tools)
        drc = distance(
                seq1 = sequtils.get_reverse_complement(seq1),
                seq2 = seq2,
                per_site = per_site,
                aligned = False,
                ignore_gaps = ignore_gaps,
                alphabet = alphabet,
                aligner_tools = aligner_tools)
        return d, drc

def distance(seq1, seq2,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle']):
    diffs, n  = get_differences(seq1 = seq1, seq2 = seq2,
            aligned = aligned,
            ignore_gaps = ignore_gaps,
            alphabet = alphabet,
            aligner_tools = aligner_tools)
    if per_site:
        return len(diffs) / float(n)
    return len(diffs)

def get_differences(seq1, seq2,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle']):
    if not alphabet:
        alphabet = DnaAlphabet()
    if not aligned:
        seq1, seq2 = align_pair(seq1, seq2, aligner_tools)
    if len(seq1) != len(seq2):
        raise AlignmentError('Sequences are not aligned')
    residue_codes = alphabet.all_residue_codes
    diffs = {}
    num_comparisons = 0
    for i in range(len(seq1)):
        if (seq1[i] == alphabet.gap) or (seq2[i] == alphabet.gap):
            if ignore_gaps:
                continue
            num_comparisons += 1
            if seq1[i].upper() != seq2[i].upper():
                diffs[i] = (seq1[i].upper(), seq2[i].upper())
            continue
        num_comparisons += 1
        if seq1[i].upper() == seq2[i].upper():
            continue
        s1 = set(residue_codes[seq1[i].upper()])
        s2 = set(residue_codes[seq2[i].upper()])
        if not s1.intersection(s2):
            diffs[i] = (seq1[i].upper(), seq2[i].upper())
    return diffs, num_comparisons

