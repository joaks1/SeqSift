import os
import sys
import unittest
import types
import itertools

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from seqsift.seqops.seqstats import *
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class ColumnFrequenciesTestCase(SeqSiftTestCase):
    def setUp(self):
        self.unaligned = [SeqRecord(Seq('ACGT'), id='1'),
                          SeqRecord(Seq('ACGTA'), id='2'),
                          SeqRecord(Seq('ACGT'), id='3'),
                          SeqRecord(Seq('ACGT'), id='4'),
                          SeqRecord(Seq('ACGT'), id='5')]
        self.simple_alignment = [
                SeqRecord(Seq('ACGT?'), id='1'),
                SeqRecord(Seq('ACGT-'), id='2'),
                SeqRecord(Seq('ACGT?'), id='3'),
                SeqRecord(Seq('ACGT-'), id='4'),
                SeqRecord(Seq('ACGT?'), id='5')]

    def test_alignment_error(self):
        self.assertRaises(AlignmentError, column_frequencies, self.unaligned)
    
    def test_unambiguous_bases(self):
        for i, char in enumerate(['A', 'C', 'G', 'T']):
            expected_freqs = [0.0] * 5
            expected_freqs[i] = 1.0
            col_freqs, seqs = column_frequencies(
                    self.simple_alignment,
                    character_list=[char])
            self.assertEqual(
                    col_freqs,
                    expected_freqs)
            self.assertSameSequences(self.simple_alignment, seqs)

    def test_gap_chars(self):
        expected_freqs = ([0.0] * 4) + [1.0]
        col_freqs, seqs = column_frequencies(
                self.simple_alignment)
        self.assertEqual(
                col_freqs,
                expected_freqs)
        self.assertSameSequences(self.simple_alignment, seqs)
        expected_freqs = ([0.0] * 4) + [3/float(5)]
        col_freqs, seqs = column_frequencies(
                self.simple_alignment,
                character_list = ['?'])
        self.assertAlmostEqual(
                col_freqs,
                expected_freqs)
        self.assertSameSequences(self.simple_alignment, seqs)
    
    def test_missing_character(self):
        expected_freqs = [0.0] * 5
        col_freqs, seqs = column_frequencies(
                self.simple_alignment,
                character_list = ['x'])
        self.assertAlmostEqual(
                col_freqs,
                expected_freqs)
        self.assertSameSequences(self.simple_alignment, seqs)

if __name__ == '__main__':
    unittest.main()

