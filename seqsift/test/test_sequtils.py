#! /usr/bin/env python

import os
import sys
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

from seqsift.seqops import sequtils
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class CopySeqMetadataTestCase(SeqSiftTestCase):
    def setUp(self):
        self.seq = SeqIO.read(
                package_paths.data_path('JF314862.gb'),
                format='gb',
                alphabet=IUPAC.ambiguous_dna)
    
    def test_empty_seq(self):
        s = sequtils.copy_seq_metadata(self.seq)
        self.assertSameMetadata(self.seq, s)
        self.assertEqual(str(s.seq), '')
        self.assertNotEqual(s, self.seq)

    def test_string(self):
        s = sequtils.copy_seq_metadata(self.seq, 'AGCT')
        self.assertSameMetadata(self.seq, s)
        self.assertEqual(str(s.seq), 'AGCT')
        self.assertNotEqual(s, self.seq)

    def test_seq(self):
        s = sequtils.copy_seq_metadata(self.seq, Seq('AGCT'))
        self.assertSameMetadata(self.seq, s)
        self.assertEqual(str(s.seq), 'AGCT')
        self.assertNotEqual(s, self.seq)

    def test_seq_record(self):
        s = sequtils.copy_seq_metadata(self.seq, SeqRecord(Seq('AGCT'), id='1'))
        self.assertSameMetadata(self.seq, s)
        self.assertEqual(str(s.seq), 'AGCT')
        self.assertNotEqual(s, self.seq)

class GetWithoutGapsTestCase(SeqSiftTestCase):
    def test_without_gaps(self):
        s1 = SeqRecord(Seq('CCATG'), id='1')
        s2 = sequtils.get_without_gaps(s1)
        self.assertSameMetadata(s1, s2)
        self.assertEqual(str(s1.seq), str(s2.seq))
        self.assertNotEqual(s1, s2)

    def test_without_gaps(self):
        s1 = SeqRecord(Seq('--CC-AT--G---'), id='1')
        s2 = sequtils.get_without_gaps(s1)
        self.assertSameMetadata(s1, s2)
        self.assertEqual(str(s2.seq), 'CCATG')
        self.assertNotEqual(s1, s2)

class GetTranslation(SeqSiftTestCase):
    def test_unambiguous_dna(self):
        s1 = SeqRecord(Seq('ATGACCAACTGA', IUPAC.unambiguous_dna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertNotEqual(s1, s2)
        self.assertEqual(str(s2.seq), 'MTN*')

    def test_unambiguous_rna(self):
        s1 = SeqRecord(Seq('AUGACCAACUGA', IUPAC.unambiguous_rna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertNotEqual(s1, s2)
        self.assertEqual(str(s2.seq), 'MTN*')

    def test_ambiguous_dna(self):
        s1 = SeqRecord(Seq('ATGATRAACTGA', IUPAC.ambiguous_dna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertNotEqual(s1, s2)
        self.assertEqual(str(s2.seq), 'MXN*')

    def test_ambiguous_rna(self):
        s1 = SeqRecord(Seq('AUGAURAACUGA', IUPAC.ambiguous_rna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertNotEqual(s1, s2)
        self.assertEqual(str(s2.seq), 'MXN*')

    def test_to_stop(self):
        s1 = SeqRecord(Seq('ATGACCAACTGA', IUPAC.unambiguous_dna), id='1')
        s2 = sequtils.get_translation(s1, to_stop = True)
        self.assertSameMetadata(s1, s2)
        self.assertNotEqual(s1, s2)
        self.assertEqual(str(s2.seq), 'MTN')

    def test_extra_base(self):
        s1 = SeqRecord(Seq('ATGACCAACTGAATAA', IUPAC.unambiguous_dna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertNotEqual(s1, s2)
        self.assertEqual(str(s2.seq), 'MTN*I')

if __name__ == '__main__':
    unittest.main()
