#! /usr/bin/env python

import os
import sys
import unittest
import types
import itertools

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from seqsift.seqops import seqmod
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class SeqModTestCase(SeqSiftTestCase):
    def setUp(self):
        self.seqs = [
                SeqRecord(Seq('ACGT'), id='1'),
                SeqRecord(Seq('ACGT'), id='2'),
                SeqRecord(Seq('ACGT'), id='3'),
                SeqRecord(Seq('ACGT'), id='4'),]
    
    def test_default(self):
        new_seqs = seqmod.seq_mod(self.seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, self.seqs)

    def test_unbalanced_strings(self):
        new_seqs = seqmod.seq_mod(self.seqs, 'A', 'GT')
        self.assertRaises(ValueError, list, new_seqs)

    def test_simple_change(self):
        new_seqs = seqmod.seq_mod(self.seqs, 'T', 'A')
        expected = [SeqRecord(str(s.seq).replace('T','A'), id=s.id) for s in self.seqs]
        self.assertSameData(new_seqs, expected)

    def test_simple_deletion(self):
        new_seqs = seqmod.seq_mod(self.seqs, del_chars='T')
        expected = [SeqRecord(str(s.seq).replace('T',''), id=s.id) for s in self.seqs]
        self.assertSameData(new_seqs, expected)

class RemoveGapsTestCase(SeqSiftTestCase):
    def setUp(self):
        self.seqs = [
                SeqRecord(Seq('ACGT-'), id='1'),
                SeqRecord(Seq('ACG-T'), id='2'),
                SeqRecord(Seq('AC-GT'), id='3'),
                SeqRecord(Seq('A-CGT'), id='4'),]

    def test_simple(self):
        new_seqs = seqmod.remove_gaps(self.seqs)
        expected = [SeqRecord(str(s.seq).replace('-',''), id=s.id) for s in self.seqs]
        self.assertSameData(new_seqs, expected)

class TranslateSeqsTestCase(SeqSiftTestCase):
    def test_dna(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ATGACCAACTGACAC', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ATGACCAACTGACCC', IUPAC.ambiguous_dna), id='3'),
                ]
        exp = [
                SeqRecord(Seq('MTN*I', IUPAC.protein), id='1'),
                SeqRecord(Seq('MTN*H', IUPAC.protein), id='2'),
                SeqRecord(Seq('MTN*P', IUPAC.protein), id='3'),
                ]
        new_seqs = seqmod.translate_seqs(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

    def test_dna_gaps(self):
        seqs = [
                SeqRecord(Seq('A-TGAC--CAACTG--AATA-', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('---ATGACCAACTGACAC---', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ATGACC-AA-C--TGACCC', IUPAC.ambiguous_dna), id='3'),
                ]
        exp = [
                SeqRecord(Seq('MTN*I', IUPAC.protein), id='1'),
                SeqRecord(Seq('MTN*H', IUPAC.protein), id='2'),
                SeqRecord(Seq('MTN*P', IUPAC.protein), id='3'),
                ]
        new_seqs = seqmod.translate_seqs(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

    def test_dna_extra_bases(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTGAATAA', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ATGACCAACTGACACTA', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ATGACCAACTGACCCAC', IUPAC.ambiguous_dna), id='3'),
                ]
        exp = [
                SeqRecord(Seq('MTN*I', IUPAC.protein), id='1'),
                SeqRecord(Seq('MTN*H', IUPAC.protein), id='2'),
                SeqRecord(Seq('MTN*P', IUPAC.protein), id='3'),
                ]
        new_seqs = seqmod.translate_seqs(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

if __name__ == '__main__':
    unittest.main()

