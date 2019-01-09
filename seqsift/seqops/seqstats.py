#! /usr/bin/env python

import sys
import os
import math
import itertools

from seqsift.align import align, align_pair
from seqsift.utils import functions, stats, dataio
from seqsift.utils.errors import AlignmentError
from seqsift.utils import alphabets
from seqsift.seqops import sequtils
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def get_duplicate_ids(seq_iter):
    dups = set()
    ids = [s.id for s in seq_iter]
    if len(ids) == len(set(ids)):
        return dups
    id_set = set()
    for i in ids:
        if i in id_set:
            dups.add(i)
        id_set.add(i)
    return sorted(list(dups))

def column_frequencies(seq_iter, character_list=['-','?']):
    seqs = dataio.BufferedIter(seq_iter)
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

def variable_columns(seq_iter):
    seqs = dataio.BufferedIter(seq_iter)
    ref_seq = None
    align_length = None
    is_variable = None
    for i, seq_record in enumerate(seqs):
        if align_length == None:
            align_length = len(seq_record)
            is_variable = [False] * align_length
            ref_seq = seq_record
            continue
        elif len(seq_record) != align_length:
            raise AlignmentError('Sequence {0} has unexpected '
                    'length {1}.'.format(seq_record.name, len(seq_record)))
        for j, character in enumerate(seq_record.seq):
            if character.lower() != ref_seq.seq[j].lower():
                is_variable[j] = True
    return is_variable, seqs

def pairwise_distance_iter(seq_iter,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle']):
    seqs = dataio.BufferedIter(seq_iter)
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

def average_number_of_pairwise_differences(seq_iter,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle']):
    sum_diffs = 0.0
    seqs = dataio.BufferedIter(seq_iter)
    n = 0
    for i, (seq1, seq2) in enumerate(itertools.combinations(seqs, 2)):
        d = distance(
                seq1 = seq1,
                seq2 = seq2,
                per_site = per_site,
                aligned = aligned,
                ignore_gaps = ignore_gaps,
                alphabet = alphabet,
                aligner_tools = aligner_tools)
        if math.isnan(d):
            continue
        sum_diffs += d
        n += 1
    if n < 1:
        return None
    return sum_diffs / n

def get_population_pair_diversity_summary(seq_iter1, seq_iter2,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle']):
    sum_diffs = 0.0
    seqs_1 = dataio.BufferedIter(seq_iter1)
    seqs_2 = dataio.BufferedIter(seq_iter2)
    pi_1 = average_number_of_pairwise_differences(seqs_1,
            per_site = per_site,
            aligned = aligned,
            ignore_gaps = ignore_gaps,
            alphabet = alphabet,
            aligner_tools = aligner_tools)
    pi_2 = average_number_of_pairwise_differences(seqs_2,
            per_site = per_site,
            aligned = aligned,
            ignore_gaps = ignore_gaps,
            alphabet = alphabet,
            aligner_tools = aligner_tools)
    assert pi_1 is not None
    assert pi_2 is not None
    n = 0
    for seq1 in seqs_1:
        for seq2 in seqs_2:
            d = distance(
                    seq1 = seq1,
                    seq2 = seq2,
                    per_site = per_site,
                    aligned = aligned,
                    ignore_gaps = ignore_gaps,
                    alphabet = alphabet,
                    aligner_tools = aligner_tools)
            if math.isnan(d):
                continue
            sum_diffs += d
            n += 1
    assert n > 0
    pi_b = sum_diffs / n
    pi_w = (pi_1 + pi_2) / 2.0
    pi_net = pi_b - pi_w
    return {
            "pi_1": pi_1,
            "pi_2": pi_2,
            "pi_within": pi_w,
            "pi_between": pi_b,
            "pi_net": pi_net,
            }

def sample_distance_iter(seq_iter,
        sample_size,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle'],
        rng = None):
    seqs = dataio.BufferedIter(seq_iter)
    seqs_to_sample = dataio.BufferedIter(seqs)
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
        rc = None
        drc = None
        if (not alphabet) or (not alphabet.has_state('M')):
            try:
                rc = sequtils.get_reverse_complement(seq1),
            except:
                pass
        if rc:
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
        if n < 1:
            return float("nan")
        return len(diffs) / float(n)
    return len(diffs)

def get_differences(seq1, seq2,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle']):
    if not alphabet:
        alphabet = alphabets.DnaAlphabet()
    if not aligned:
        seq1, seq2 = align_pair(seq1, seq2, aligner_tools)
    if len(seq1) != len(seq2):
        raise AlignmentError('Sequences are not aligned')
    residue_codes = alphabet.all_residue_codes
    diffs = {}
    num_comparisons = 0
    for i in range(len(seq1)):
        if (seq1[i] == alphabet.missing) or (seq2[i] == alphabet.missing):
            continue
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

def get_seq_summary(seqs):
    return stats.SampleSummarizer((len(s) for s in seqs))

def get_seq_summaries_from_files(paths, format = None, data_type = 'dna',
        ambiguities = True, global_key = 'global'):
    seq_iters = (dataio.SeqFileIter(p,
            format = format,
            data_type = data_type,
            ambiguities = ambiguities
            ) for p in paths)
    sums = {}
    sums[global_key] = stats.SampleSummarizer()
    for i, seq_iter in enumerate(seq_iters):
        key = seq_iter.name
        if key == global_key:
            raise Exception('{0} is the global key for the seq summary '
                    'dict.\nSpecify alternative global key'.format(key))
        sums[key] = stats.SampleSummarizer()
        for seq in seq_iter:
            sums[key].add_sample(len(seq))
            sums[global_key].add_sample(len(seq))
    return sums

