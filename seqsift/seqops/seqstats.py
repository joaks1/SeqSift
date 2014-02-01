#! /usr/bin/env python

import sys
import os
import itertools
import datetime

from seqsift.align import align, align_pair
from seqsift.utils import functions, stats
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

def summarize_distances(seq_iter,
        sample_size = 0,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        do_full_alignment = False,
        full_alignment_out_path = None,
        aligner_tools = ['mafft', 'muscle'],
        full_aligner_tools = None,
        rng = None,
        log_frequency = 0):
    if ((not aligned) and (do_full_alignment)):
        if not full_aligner_tools:
            full_aligner_tools = aligner_tools
        seq_iter = align(seq_iter,
                tools = full_aligner_tools,
                out_path = full_alignment_out_path)
        aligned = True
    if sample_size > 0:
        distance_iter = sample_distance_iter(
                seq_iter = seq_iter,
                sample_size = sample_size,
                aligned = aligned,
                ignore_gaps = ignore_gaps,
                per_site = per_site,
                alphabet = alphabet,
                aligner_tools = aligner_tools,
                rng = rng)
    else:
        distance_iter = pairwise_distance_iter(
                seq_iter = seq_iter,
                aligned = aligned,
                ignore_gaps = ignore_gaps,
                per_site = per_site,
                alphabet = alphabet,
                aligner_tools = aligner_tools)
    distances = {}
    rev_comp_errors = []
    for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
        if (log_frequency > 0) and (((i + 1) % log_frequency) == 0):
            _LOG.info('{0}: Calulating distance for comparison {1}...'.format(
                    datetime.datetime.now(),
                    (i + 1)))
        if drc < d:
            _LOG.warning('reverse complement of {0} is more similar to '
                    '{1}'.format(seq1.id, seq2.id))
            rev_comp_errors.append((seq1.id, seq2.id, d, drc))
        if sample_size > 0:
            if not distances.has_key(seq1.id):
                distances[seq1.id] = stats.SampleSummarizer(samples = [d])
                continue
            distances[seq1.id].add_sample(d)
        else:
            if not distances.has_key(seq1.id):
                distances[seq1.id] = stats.SampleSummarizer(samples = [d])
            else:
                distances[seq1.id].add_sample(d)
            if not distances.has_key(seq2.id):
                distances[seq2.id] = stats.SampleSummarizer(samples = [d])
            else:
                distances[seq2.id].add_sample(d)
    return distances, rev_comp_errors

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
    diffs, l  = get_differences(seq1 = seq1, seq2 = seq2,
            aligned = aligned,
            ignore_gaps = ignore_gaps,
            alphabet = alphabet,
            aligner_tools = aligner_tools)
    if per_site:
        return len(diffs) / float(l)
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

