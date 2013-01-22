#! /usr/bin/env python

import sys
import os
from string import maketrans

from Bio.SeqRecord import SeqRecord

from seqsift.seqops.sequtils import copy_seq_metadata
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def translate_seqs(seq_iter,
        gap_characters=['-'],
        table='Standard',
        stop_symbol='*',
        to_stop=False,
        cds=False):
    '''
    Translates DNA or RNA sequences into amino acid sequences.

    This functions calls the Seq.translate method of biopython and returns a
    generator function of resulting amino acid sequences.

    The arguments are passed direclty to Seq.translate. Below is biopython's
    description of the arguments.

    Arguments:
     - table - Which codon table to use?  This can be either a name
               (string), an NCBI identifier (integer), or a CodonTable
               object (useful for non-standard genetic codes).  This
               defaults to the "Standard" table.
     - stop_symbol - Single character string, what to use for terminators.
                     This defaults to the asterisk, "*".
     - to_stop - Boolean, defaults to False meaning do a full translation
                 continuing on past any stop codons (translated as the
                 specified stop_symbol).  If True, translation is
                 terminated at the first in frame stop codon (and the
                 stop_symbol is not appended to the returned protein
                 sequence).
     - cds - Boolean, indicates this is a complete CDS.  If True,
             this checks the sequence starts with a valid alternative start
             codon (which will be translated as methionine, M), that the
             sequence length is a multiple of three, and that there is a
             single in frame stop codon at the end (this will be excluded
             from the protein sequence, regardless of the to_stop option).
             If these tests fail, an exception is raised.
    '''
    for s in remove_gaps(seq_iter, gap_characters=gap_characters):
        yield SeqRecord(
                seq = s.seq.translate(table=table, stop_symbol=stop_symbol,
                        to_stop=to_stop, cds=cds),
                id = s.id,
                name = s.name,
                description = s.description)

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
