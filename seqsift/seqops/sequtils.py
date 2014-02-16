#! /usr/bin/env python

import sys
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def copy_seq_metadata(seq_record, new_seq=''):
    if isinstance(new_seq, Seq):
        s = new_seq
    elif hasattr(new_seq, 'seq'):
        s = new_seq.seq
    else:
        s = Seq(new_seq, alphabet=seq_record.seq.alphabet)
    return SeqRecord(
            seq = s,
            id = seq_record.id,
            name = seq_record.name,
            description = seq_record.description,
            letter_annotations = seq_record.letter_annotations,
            annotations = seq_record.annotations,
            features = seq_record.features,
            dbxrefs = seq_record.dbxrefs)

def sequences_are_equal(seq_record1, seq_record2):
    if seq_record1.id != seq_record2.id:
        return False
    if seq_record1.name != seq_record2.name:
        return False
    if seq_record1.description != seq_record2.description:
        return False
    if str(seq_record1.seq) != str(seq_record2.seq):
        return False
    return True

def get_reverse_complement(seq_record):
    return copy_seq_metadata(seq_record,
            new_seq = str(seq_record.seq.reverse_complement()))

def get_without_gaps(seq_record):
    return copy_seq_metadata(seq_record,
            new_seq = ''.join([x for x in str(seq_record.seq) if x != '-']))

def get_translation(seq_record, **kwargs):
    if not hasattr(seq_record.seq, 'translate'):
        raise Exception('seq record {0!r} does not have a translate '
                'method'.format(seq_record))
    s = seq_record.seq
    while len(s) % 3 != 0:
        s = Seq(str(s)[:-1], alphabet = s.alphabet)
    return copy_seq_metadata(seq_record,
            new_seq = s.translate(**kwargs))

def get_longest_reading_frames(seq_record, table = 1,
        allow_partial = True,
        require_start_after_stop = True):
    lrf = []
    for i in range(3):
        p = get_translation(
                copy_seq_metadata(
                        seq_record, new_seq = seq_record.seq[i:]),
                table = table,
                to_stop = False)
        fragments = str(p.seq).split('*')
        if (not allow_partial) and (not str(p.seq).endswith('*')):
            fragments.pop(-1)
        if len(fragments) < 1:
            continue
        stop_index = i
        start_offset = 0
        for j, f in enumerate(fragments):
            start_index = stop_index
            stop_index += (3 * len(f)) + 3
            if len(f) < 1:
                continue
            if (not allow_partial) or (require_start_after_stop and (j > 0)):
                start_offset = 3 * f.find('M')
                if start_offset < 0:
                    continue
            if stop_index > (start_index + start_offset):
                rf = copy_seq_metadata(seq_record,
                            new_seq = seq_record.seq[
                                    (start_index + start_offset): stop_index])
                if len(rf.seq) < 1:
                    continue
                if len(lrf) < 1:
                    lrf = [rf]
                    continue
                if len(rf.seq) == len(lrf[0].seq):
                    lrf.append(rf)
                    continue
                if len(rf.seq) > len(lrf[0].seq):
                    lrf = [rf]
                    continue
    return lrf

