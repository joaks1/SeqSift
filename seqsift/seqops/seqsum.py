#! /usr/bin/env python

import sys
import os
import datetime

from seqsift import align
from seqsift.utils import stats
from seqsift.seqops import seqstats, seqmod
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def summarize_longest_read_lengths(seq_iter,
        gap_characters=['-'],
        table = 1,
        allow_partial = True,
        require_start_after_stop = True):
    lengths = []
    for seq in seqmod.longest_reading_frames(seq_iter,
            gap_characters = gap_characters,
            table = table,
            allow_partial = allow_partial,
            require_start_after_stop = require_start_after_stop):
        lengths.append((len(seq.seq), seq.id))
    return sorted(lengths)

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
        seq_iter = align.align(seq_iter,
                tools = full_aligner_tools,
                out_path = full_alignment_out_path)
        aligned = True
    if sample_size > 0:
        distance_iter = seqstats.sample_distance_iter(
                seq_iter = seq_iter,
                sample_size = sample_size,
                aligned = aligned,
                ignore_gaps = ignore_gaps,
                per_site = per_site,
                alphabet = alphabet,
                aligner_tools = aligner_tools,
                rng = rng)
    else:
        distance_iter = seqstats.pairwise_distance_iter(
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
                    '{1} ({2:.5f} vs {3:.5f})'.format(seq1.id, seq2.id, drc, d))
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

