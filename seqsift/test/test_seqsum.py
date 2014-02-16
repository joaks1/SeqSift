#! /usr/bin/env python

import os
import sys
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqsift.seqops import seqsum
from seqsift.utils import dataio
from seqsift.utils import functions, errors
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class SummarizeLongestReadLengthsTestCase(unittest.TestCase):
    def test_default(self):
        seqs = [
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACAAACTCCTCCTGAATA'), id='1'),
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACAAACTCCTCCTCCTGAATA'), id='2'),
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACAAACTGAATA'), id='3'),
                ]
        exp = [
                (12, 35, '3'),
                (18, 41, '1'),
                (21, 44, '2'),
                ]
        lengths = seqsum.summarize_longest_read_lengths(seqs,
                gap_characters=['-'],
                table = 1,
                allow_partial = True,
                require_start_after_stop = True)
        self.assertEqual(lengths, exp)

    def test_no_partial(self):
        seqs = [
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACAAACTCCTCCTGAATA'), id='1'),
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACAAACTCCTCCTCCTGAATA'), id='2'),
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACAAACTGAATA'), id='3'),
                ]
        exp = [
                (12, 0, '3'),
                (18, 0, '1'),
                (21, 0, '2'),
                ]
        lengths = seqsum.summarize_longest_read_lengths(seqs,
                gap_characters=['-'],
                table = 1,
                allow_partial = False,
                require_start_after_stop = True)
        self.assertEqual(lengths, exp)

class SummarizeDistancesTestCase(unittest.TestCase):
    def setUp(self):
        self.seqs = [
                SeqRecord(Seq('CCA--CGTAA'), id='1'),
                SeqRecord(Seq('CCG--CGTAA'), id='2'),
                SeqRecord(Seq('CCA--TATAA'), id='3')]
        self.expected_means = {'1': 1.5,
                               '2': 2.0,
                               '3': 2.5}
        self.expected_maxs = {'1': 2,
                               '2': 3,
                               '3': 3}
        self.rc_path = package_paths.data_path('primates-rev-comp-error.fasta')
        self.gappy_path = package_paths.data_path('melittobia-its1.fasta')
        self.rc_gappy_path = package_paths.data_path(
                'melittobia-its1-rev-comp-error.fasta')


    def test_aligned(self):
        d, e = seqsum.summarize_distances(self.seqs,
                sample_size = 0,
                per_site = False,
                aligned = True,
                ignore_gaps = True)
        self.assertEqual(e, [])
        self.assertEqual(sorted(d.keys()), sorted(self.expected_means.keys()))
        for k in d.iterkeys():
            self.assertEqual(d[k].maximum, self.expected_maxs[k])
            self.assertAlmostEqual(d[k].mean, self.expected_means[k])

    def test_unaligned(self):
        d, e = seqsum.summarize_distances(self.seqs,
                sample_size = 0,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = None)
        self.assertEqual(e, [])
        self.assertEqual(sorted(d.keys()), sorted(self.expected_means.keys()))
        for k in d.iterkeys():
            self.assertEqual(d[k].maximum, self.expected_maxs[k])
            self.assertAlmostEqual(d[k].mean, self.expected_means[k])

    def test_unaligned_mafft(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        d, e = seqsum.summarize_distances(self.seqs,
                sample_size = 0,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = ['mafft'])
        self.assertEqual(e, [])
        self.assertEqual(sorted(d.keys()), sorted(self.expected_means.keys()))
        for k in d.iterkeys():
            self.assertEqual(d[k].maximum, self.expected_maxs[k])
            self.assertAlmostEqual(d[k].mean, self.expected_means[k])

    def test_unaligned_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        d, e = seqsum.summarize_distances(self.seqs,
                sample_size = 0,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = ['muscle'])
        self.assertEqual(e, [])
        self.assertEqual(sorted(d.keys()), sorted(self.expected_means.keys()))
        for k in d.iterkeys():
            self.assertEqual(d[k].maximum, self.expected_maxs[k])
            self.assertAlmostEqual(d[k].mean, self.expected_means[k])

    def test_full_alignment_mafft(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        d, e = seqsum.summarize_distances(self.seqs,
                sample_size = 0,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = True,
                aligner_tools = ['mafft'])
        self.assertEqual(e, [])
        self.assertEqual(sorted(d.keys()), sorted(self.expected_means.keys()))
        for k in d.iterkeys():
            self.assertEqual(d[k].maximum, self.expected_maxs[k])
            self.assertAlmostEqual(d[k].mean, self.expected_means[k])

    def test_full_alignment_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        d, e = seqsum.summarize_distances(self.seqs,
                sample_size = 0,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = True,
                aligner_tools = ['muscle'])
        self.assertEqual(e, [])
        self.assertEqual(sorted(d.keys()), sorted(self.expected_means.keys()))
        for k in d.iterkeys():
            self.assertEqual(d[k].maximum, self.expected_maxs[k])
            self.assertAlmostEqual(d[k].mean, self.expected_means[k])

    def test_full_alignment_error(self):
        self.assertRaises(errors.ExternalToolNotFoundError,
                seqsum.summarize_distances,
                self.seqs,
                0,
                False,
                False,
                True,
                None,
                True,
                None,
                None,
                None)

    def test_rev_comp_error_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqsum.summarize_distances(self.rc_seqs,
                sample_size = 0,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = ['muscle'],
                full_aligner_tools = None)
        self.assertEqual(len(e), 11)
        for rce in e:
            self.assertTrue('Homo_sapiens' in rce)
        self.assertEqual(len(d), 12)

    def test_rev_comp_error_muscle_sample(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqsum.summarize_distances(self.rc_seqs,
                sample_size = 5,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = ['muscle'],
                full_aligner_tools = None)
        self.assertTrue(len(e) >= 5)
        for rce in e:
            self.assertTrue('Homo_sapiens' in rce)
        self.assertEqual(len(d), 12)

    def test_rev_comp_error_muscle_full(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqsum.summarize_distances(self.rc_seqs,
                sample_size = 0,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = True,
                aligner_tools = ['muscle'],
                full_aligner_tools = ['muscle'])
        self.assertEqual(len(e), 11)
        for rce in e:
            self.assertTrue('Homo_sapiens' in rce)
        self.assertEqual(len(d), 12)

    def test_rev_comp_error_mafft(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqsum.summarize_distances(self.rc_seqs,
                sample_size = 0,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = ['mafft'],
                full_aligner_tools = None)
        self.assertEqual(len(e), 11)
        for rce in e:
            self.assertTrue('Homo_sapiens' in rce)
        self.assertEqual(len(d), 12)

    def test_rev_comp_error_mafft_sample(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqsum.summarize_distances(self.rc_seqs,
                sample_size = 5,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = ['mafft'],
                full_aligner_tools = None)
        self.assertTrue(len(e) >= 5)
        for rce in e:
            self.assertTrue('Homo_sapiens' in rce)
        self.assertEqual(len(d), 12)

    def test_rev_comp_error_mafft_full(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqsum.summarize_distances(self.rc_seqs,
                sample_size = 0,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = True,
                aligner_tools = ['mafft'],
                full_aligner_tools = ['mafft'])
        self.assertEqual(len(e), 11)
        for rce in e:
            self.assertTrue('Homo_sapiens' in rce)
        self.assertEqual(len(d), 12)

    def test_rev_comp_gappy_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.gappy_path)
        d, e = seqsum.summarize_distances(self.rc_seqs,
                sample_size = 5,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = ['muscle'],
                full_aligner_tools = None)
        self.assertTrue(len(e) < 1)
        self.assertEqual(len(d), 31)

    def test_rev_comp_error_gappy_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_gappy_path)
        d, e = seqsum.summarize_distances(self.rc_seqs,
                sample_size = 5,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = ['muscle'],
                full_aligner_tools = None)
        self.assertTrue(len(e) >= 5)
        for rce in e:
            self.assertTrue('JF924943_Dibrachys_pelos' in rce)
        self.assertEqual(len(d), 31)

if __name__ == '__main__':
    unittest.main()

