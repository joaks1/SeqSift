#! /usr/bin/env python

import sys
import os
import datetime
from string import maketrans

from Bio.SeqRecord import SeqRecord

from seqsift.seqops import sequtils
from seqsift.seqops import seqstats
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def longest_reading_frames(seq_iter,
        gap_characters=['-'],
        table = 1,
        allow_partial = True,
        require_start_after_stop = True):
    for i, s in enumerate(remove_gaps(seq_iter, gap_characters=gap_characters)):
        lrf = sequtils.get_longest_reading_frames(s,
                table = table,
                allow_partial = allow_partial,
                require_start_after_stop = require_start_after_stop)
        if lrf:
            yield lrf[0]
        else:
            yield sequtils.copy_seq_metadata(s, '')


def translate_longest_reading_frames(seq_iter,
        gap_characters=['-'],
        table = 1,
        allow_partial = True,
        require_start_after_stop = True):
    frames = longest_reading_frames(seq_iter,
            gap_characters = gap_characters,
            table = table,
            allow_partial = allow_partial,
            require_start_after_stop = require_start_after_stop)
    return translate_seqs(frames, table = table)

def reverse_complement(seq_iter):
    for s in seq_iter:
        yield sequtils.get_reverse_complement(s)

def reverse_complement_to_first_seq(seq_iter,
        per_site = True,
        aligned = False,
        ignore_gaps = True,
        alphabet = None,
        aligner_tools = ['mafft', 'muscle'],
        log_frequency = 0):
    seq1 = None
    for i, seq2 in enumerate(seq_iter):
        if i == 0:
            seq1 = seq2
            yield seq2
            continue
        if (log_frequency > 0) and (((i + 1) % log_frequency) == 0):
            _LOG.info('{0}: Checking reverse complement of seq {1}...'.format(
                    datetime.datetime.now(),
                    (i + 1)))
        d, drc = seqstats.get_distances(
                seq1 = seq1,
                seq2 = seq2,
                per_site = per_site,
                aligned = aligned,
                ignore_gaps = ignore_gaps,
                alphabet = alphabet,
                aligner_tools = aligner_tools)
        _LOG.debug('{0}: distance {1}, rev comp distance {2}'.format(
                seq2.id, d, drc))
        if drc < d:
            _LOG.warning('Reverse complementing sequence {0!r} (length {1})\n\t'
                    'rev comp distance ({2}) < current distance '
                    '({3})'.format(seq2.id, len(seq2.seq), drc, d))
            yield sequtils.get_reverse_complement(seq2)
            continue
        yield seq2

def reverse_complement_to_longest_reading_frame(seq_iter,
        gap_characters=['-'],
        table = 1,
        allow_partial = True,
        require_start_after_stop = True,
        log_frequency = 0):
    for i, s in enumerate(remove_gaps(seq_iter, gap_characters=gap_characters)):
        if (log_frequency > 0) and (((i + 1) % log_frequency) == 0):
            _LOG.info('{0}: Checking reverse complement of seq {1}...'.format(
                    datetime.datetime.now(),
                    (i + 1)))
        rc = sequtils.get_reverse_complement(s)
        p1 = sequtils.get_longest_reading_frames(seq_record = s,
                table = table,
                allow_partial = allow_partial,
                require_start_after_stop = require_start_after_stop)
        p2 = sequtils.get_longest_reading_frames(seq_record = rc,
                table = table,
                allow_partial = allow_partial,
                require_start_after_stop = require_start_after_stop)
        _LOG.debug('{0}: read length {1}, rev comp read length {2}'.format(
                s.id, len(p1[0].seq), len(p2[0].seq)))
        if len(p2) == 0:
            yield s
        elif len(p1) == 0:
            _LOG.warning('Reverse complementing sequence {0!r}'.format(rc.id))
            yield rc
        elif len(p2[0].seq) > len(p1[0].seq):
            _LOG.warning('Reverse complementing sequence {0!r}'.format(rc.id))
            yield rc
        else:
            yield s

def translate_seqs(seq_iter,
        gap_characters=['-'],
        **kwargs):
    '''
    Translates DNA or RNA sequences into amino acid sequences.

    This function returns a generator that yields a copy of each input sequence
    translated into an amino acid sequence (all `gap_characters` are removed
    from the sequence prior to translation). `kwargs` are keyword arguments
    that are passed to the Seq.translate method of biopython to translate each
    sequence. Below is biopython's description of valid keyword arguments:

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
        yield sequtils.get_translation(seq_record = s, **kwargs)

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
        yield sequtils.copy_seq_metadata(seq,
                new_seq=str(seq.seq).translate(table, del_chars))

def dice(seq_iter,
        slices_to_keep):
    for seq in seq_iter:
        yield sequtils.copy_seq_metadata(seq,
                new_seq = ''.join((str(seq.seq[l:r]) for l,r in slices_to_keep)))

