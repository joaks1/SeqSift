#! /usr/bin/env python

import sys
import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def copy_seq_metadata(seq_record, new_seq=''):
    return SeqRecord(
            seq = Seq(new_seq, alphabet=seq_record.seq.alphabet),
            id = seq_record.id,
            name = seq_record.name,
            description = seq_record.description,
            letter_annotations = seq_record.letter_annotations,
            annotations = seq_record.annotations,
            features = seq_record.features,
            dbxrefs = seq_record.dbxrefs)

def get_reverse_complement(seq_record):
    return copy_seq_metadata(seq_record,
            new_seq = str(seq_record.seq.reverse_complement()))


def get_without_gaps(seq_record):
    return copy_seq_metadata(seq_record,
            new_seq = ''.join([x for x in str(seq_record.seq) if x != '-']))

