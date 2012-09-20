#! /usr/bin/env python

import os
import sys
import unittest
import types
import itertools

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from seqsift.seqops.seqfilter import *
from seqsift.utils.errors import AlignmentError
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class LengthFilterTestCase(unittest.TestCase):
    def setUp(self):
        self.seqs = [
                SeqRecord(Seq('ACGT'), id='4'),
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

class ColumnFilterTestCase(SeqSiftTestCase):
    def setUp(self):
        self.unaligned = [
                SeqRecord(Seq('ACGT'), id='1'),
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
        f = column_filter(self.unaligned)
        self.assertRaises(AlignmentError, f.next)

    def test_remove_gap(self):
        expected = [
                SeqRecord(Seq('ACGT'), id='1'),
                SeqRecord(Seq('ACGT'), id='2'),
                SeqRecord(Seq('ACGT'), id='3'),
                SeqRecord(Seq('ACGT'), id='4'),
                SeqRecord(Seq('ACGT'), id='5')]

        seqs = list(column_filter(self.simple_alignment,
                character_list = ['?', '-'],
                max_frequency = 1.0))
        self.assertSameData(seqs, expected)

        seqs = list(column_filter(self.simple_alignment,
                character_list = ['-'],
                max_frequency = 2/float(5)))
        self.assertSameData(seqs, expected)

        seqs = list(column_filter(self.simple_alignment,
                character_list = ['?'],
                max_frequency = 3/float(5)))
        self.assertSameData(seqs, expected)

        seqs = list(column_filter(self.simple_alignment,
                character_list = ['-'],
                max_frequency = 3/float(5)))
        self.assertSameData(seqs, self.simple_alignment)

        seqs = list(column_filter(self.simple_alignment,
                character_list = ['?'],
                max_frequency = 3.1/float(5)))
        self.assertSameData(seqs, self.simple_alignment)

class RowFilterTestCase(SeqSiftTestCase):
    def setUp(self): 
        self.unaligned = [
                SeqRecord(Seq('ACGT'), id='1'),
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
        
    def test_clean(self):
        seqs = list(row_filter(self.unaligned,
                character_list=['?', '-'],
                max_frequency=0.5))
        self.assertSameData(seqs, self.unaligned)
        seqs = list(row_filter(self.simple_alignment,
                character_list=['?', '-'],
                max_frequency=0.5))
        self.assertSameData(seqs, self.simple_alignment)

    def test_simple(self):
        seqs = list(row_filter(self.unaligned,
            character_list=['a'],
            max_frequency = 2/float(5)))
        self.assertSameData(seqs, self.unaligned[:1] + self.unaligned[2:])

        seqs = list(row_filter(self.simple_alignment,
            character_list=['-'],
            max_frequency = 1/float(5)))
        self.assertSameData(seqs, self.simple_alignment[::2])
        
    def test_no_filtrate(self):
        seqs = list(row_filter(self.simple_alignment,
            character_list=['?','-'],
            max_frequency = 1/float(5)))
        self.assertEqual(seqs, [])

if __name__ == '__main__':
    unittest.main()
