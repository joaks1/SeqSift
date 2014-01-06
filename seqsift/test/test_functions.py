#! /usr/bin/env python

import os
import sys
import unittest

from seqsift.utils import functions
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class TestSampleIter(unittest.TestCase):
    def test_sample_size_error(self):
        i = (x for x in range(10))
        self.assertRaises(ValueError, functions.sample_iter, i, 11)

    def test_sample_iter(self):
        counts = [0 for i in range(10)]
        for i in range(100000):
            iterable = (x for x in range(10))
            s = functions.sample_iter(iterable, 5)
            for j in s:
                counts[j] += 1
        total = sum(counts)
        self.assertEqual(total, 5 * 100000)
        freqs = [i/float(total) for i in counts]
        for f in freqs:
            self.assertAlmostEqual(f, 1/float(10), 2)

    def test_sample_iter_exclude(self):
        counts = [0 for i in range(10)]
        for i in range(100000):
            iterable = (x for x in range(10))
            s = functions.sample_iter(iterable, 5, exclude = [9])
            for j in s:
                counts[j] += 1
        total = sum(counts)
        self.assertEqual(total, 5 * 100000)
        self.assertEqual(counts[9], 0)
        freqs = [i/float(total) for i in counts if i > 0]
        self.assertEqual(len(freqs), 9)
        for f in freqs:
            self.assertAlmostEqual(f, 1/float(9), 2)

        counts = [0 for i in range(10)]
        for i in range(100000):
            iterable = (x for x in range(10))
            s = functions.sample_iter(iterable, 5, exclude = [8, 9])
            for j in s:
                counts[j] += 1
        total = sum(counts)
        self.assertEqual(total, 5 * 100000)
        self.assertEqual(counts[9], 0)
        self.assertEqual(counts[8], 0)
        freqs = [i/float(total) for i in counts if i > 0]
        self.assertEqual(len(freqs), 8)
        for f in freqs:
            self.assertAlmostEqual(f, 1/float(8), 2)

class GetNewPathTestCase(SeqSiftTestCase):
    def setUp(self):
        self.unique = os.path.join(self.temp_fs.base_dir, 
                self.temp_fs.token_id + '-blah.txt')

    def test_unique(self):
        p = functions.get_new_path(self.unique)
        self.assertEqual(p, self.unique)
        self.assertTrue(os.path.exists(p))

    def test_exists(self):
        f = open(self.unique, 'w')
        f.close()
        p = functions.get_new_path(self.unique)
        new_p = self.unique + '-0'
        self.assertEqual(p, new_p)
        self.assertTrue(os.path.exists(p))

if __name__ == '__main__':
    unittest.main()
