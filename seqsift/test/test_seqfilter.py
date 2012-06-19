import os
import sys
import unittest
import types
import itertools

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from seqsift.seqops.seqfilter import *
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class LengthFilterTestCase(unittest.TestCase):
    def setUp(self):
        self.seqs = [SeqRecord(Seq('ACGT'), id='4'),
                     SeqRecord(Seq('ACGTA'), id='5'),
                     SeqRecord(Seq('ACGTAC'), id='6'),
                     SeqRecord(Seq('ACGTACG'), id='7'),
                     SeqRecord(Seq('ACGTACGT'), id='8')]
    
    def test_limit_error(self):
        f = length_filter([], 6, 5)
        self.assertRaises(ValueError, f.next)

    def test_default(self):
        filter_iter = length_filter(self.seqs)
        self.assertIsInstance(filter_iter, types.GeneratorType)
        filtrate = list(filter_iter)
        self.assertEqual(len(filtrate), len(self.seqs))
        self.assertEqual(self.seqs, filtrate)

    def test_no_filter_limit(self):
        filtrate = list(length_filter(self.seqs, 4, 8))
        self.assertEqual(len(filtrate), len(self.seqs))
        self.assertEqual(self.seqs, filtrate)

    def test_no_max(self):
        for i in range(4, 9):
            filtrate = list(length_filter(self.seqs, i))
            self.assertEqual(len(filtrate), 9-i)

    def test_no_min(self):
        for i in range(4, 9):
            filtrate = list(length_filter(self.seqs, max_length=i))
            self.assertEqual(len(filtrate), i-3)

    def test_equal_min_max(self):
        filtrate = list(length_filter(self.seqs, 6, 6))
        self.assertEqual(len(filtrate), 1)
        self.assertEqual(filtrate[0], self.seqs[2])
    
    def test_empty_iter(self):
        filtrate = list(length_filter([]))
        self.assertEqual([], filtrate)

    def test_no_filtrate(self):
        filtrate = list(length_filter(self.seqs, 9, 10))
        self.assertEqual([], filtrate)
        filtrate = list(length_filter(self.seqs, 0, 3))
        self.assertEqual([], filtrate)

if __name__ == '__main__':
    unittest.main()
