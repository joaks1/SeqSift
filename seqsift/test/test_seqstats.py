#! /usr/bin/env python

import os
import sys
import unittest
import types
import itertools

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from seqsift.seqops import seqstats
from seqsift.utils import dataio
from seqsift.utils import functions, errors
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class TestSummarizeDistancesTestCase(unittest.TestCase):
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
        d, e = seqstats.summarize_distances(self.seqs,
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
        d, e = seqstats.summarize_distances(self.seqs,
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
        d, e = seqstats.summarize_distances(self.seqs,
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
        d, e = seqstats.summarize_distances(self.seqs,
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
        d, e = seqstats.summarize_distances(self.seqs,
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
        d, e = seqstats.summarize_distances(self.seqs,
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
                seqstats.summarize_distances,
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
        d, e = seqstats.summarize_distances(self.rc_seqs,
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

    def test_rev_comp_error_muscle_sample(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqstats.summarize_distances(self.rc_seqs,
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

    def test_rev_comp_error_muscle_full(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqstats.summarize_distances(self.rc_seqs,
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

    def test_rev_comp_error_mafft(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqstats.summarize_distances(self.rc_seqs,
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

    def test_rev_comp_error_mafft_sample(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqstats.summarize_distances(self.rc_seqs,
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

    def test_rev_comp_error_mafft_full(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_path)
        d, e = seqstats.summarize_distances(self.rc_seqs,
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

    def test_rev_comp_gappy_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.gappy_path)
        d, e = seqstats.summarize_distances(self.rc_seqs,
                sample_size = 5,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                do_full_alignment = False,
                aligner_tools = ['muscle'],
                full_aligner_tools = None)
        self.assertTrue(len(e) < 1)

    def test_rev_comp_error_gappy_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        self.rc_seqs = dataio.get_buffered_seq_iter(self.rc_gappy_path)
        d, e = seqstats.summarize_distances(self.rc_seqs,
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

class PairwiseDistanceIterTestCase(unittest.TestCase):
    def setUp(self):
        self.seqs = [
                SeqRecord(Seq('A--CGT'), id='1'),
                SeqRecord(Seq('G--CGT'), id='2'),
                SeqRecord(Seq('A--TAT'), id='3')]
        self.expected = {}
        self.expected['1'] = {'2': 1, '3': 2}
        self.expected['2'] = {'1': 1, '3': 3}
        self.expected['3'] = {'1': 2, '2': 3}

    def test_aligned(self):
        distance_iter = seqstats.pairwise_distance_iter(
                seq_iter = self.seqs,
                per_site = False,
                aligned = True,
                ignore_gaps = True)
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertEqual(
                    self.expected[seq1.id][seq2.id],
                    d)
        self.assertEqual(i, 2)

        distance_iter = seqstats.pairwise_distance_iter(
                seq_iter = self.seqs,
                per_site = True,
                aligned = True,
                ignore_gaps = True)
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertAlmostEqual(
                    self.expected[seq1.id][seq2.id] / 6.0,
                    d)
        self.assertEqual(i, 2)

    def test_unaligned(self):
        distance_iter = seqstats.pairwise_distance_iter(
                seq_iter = self.seqs,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = None)
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertEqual(
                    self.expected[seq1.id][seq2.id],
                    d)
        self.assertEqual(i, 2)

        distance_iter = seqstats.pairwise_distance_iter(
                seq_iter = self.seqs,
                per_site = True,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = None)
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertAlmostEqual(
                    self.expected[seq1.id][seq2.id] / 4.0,
                    d)
        self.assertEqual(i, 2)

    def test_unaligned_mafft(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        distance_iter = seqstats.pairwise_distance_iter(
                seq_iter = self.seqs,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = ['mafft'])
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertEqual(
                    self.expected[seq1.id][seq2.id],
                    d)
        self.assertEqual(i, 2)

        distance_iter = seqstats.pairwise_distance_iter(
                seq_iter = self.seqs,
                per_site = True,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = ['mafft'])
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertAlmostEqual(
                    self.expected[seq1.id][seq2.id] / 4.0,
                    d)
        self.assertEqual(i, 2)

    def test_unaligned_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        distance_iter = seqstats.pairwise_distance_iter(
                seq_iter = self.seqs,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = ['muscle'])
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertEqual(
                    self.expected[seq1.id][seq2.id],
                    d)
        self.assertEqual(i, 2)

        distance_iter = seqstats.pairwise_distance_iter(
                seq_iter = self.seqs,
                per_site = True,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = ['muscle'])
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertAlmostEqual(
                    self.expected[seq1.id][seq2.id] / 4.0,
                    d)
        self.assertEqual(i, 2)

class SampleDistanceIterTestCase(unittest.TestCase):
    def setUp(self):
        self.seqs = [
                SeqRecord(Seq('A--CGT'), id='1'),
                SeqRecord(Seq('G--CGT'), id='2'),
                SeqRecord(Seq('A--TAT'), id='3')]
        self.expected = {}
        self.expected['1'] = {'2': 1, '3': 2}
        self.expected['2'] = {'1': 1, '3': 3}
        self.expected['3'] = {'1': 2, '2': 3}

    def test_aligned(self):
        distance_iter = seqstats.sample_distance_iter(
                seq_iter = self.seqs,
                sample_size = 2,
                per_site = False,
                aligned = True,
                ignore_gaps = True)
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertEqual(
                    self.expected[seq1.id][seq2.id],
                    d)
        self.assertEqual(i, 5)

        distance_iter = seqstats.sample_distance_iter(
                seq_iter = self.seqs,
                sample_size = 2,
                per_site = True,
                aligned = True,
                ignore_gaps = True)
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertAlmostEqual(
                    self.expected[seq1.id][seq2.id] / 6.0,
                    d)
        self.assertEqual(i, 5)

    def test_unaligned(self):
        distance_iter = seqstats.sample_distance_iter(
                seq_iter = self.seqs,
                sample_size = 2,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = None)
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertEqual(
                    self.expected[seq1.id][seq2.id],
                    d)
        self.assertEqual(i, 5)

        distance_iter = seqstats.sample_distance_iter(
                seq_iter = self.seqs,
                sample_size = 2,
                per_site = True,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = None)
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertAlmostEqual(
                    self.expected[seq1.id][seq2.id] / 4.0,
                    d)
        self.assertEqual(i, 5)

    def test_unaligned_mafft(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        distance_iter = seqstats.sample_distance_iter(
                seq_iter = self.seqs,
                sample_size = 2,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = ['mafft'])
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertEqual(
                    self.expected[seq1.id][seq2.id],
                    d)
        self.assertEqual(i, 5)

        distance_iter = seqstats.sample_distance_iter(
                seq_iter = self.seqs,
                sample_size = 2,
                per_site = True,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = ['mafft'])
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertAlmostEqual(
                    self.expected[seq1.id][seq2.id] / 4.0,
                    d)
        self.assertEqual(i, 5)

    def test_unaligned_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        distance_iter = seqstats.sample_distance_iter(
                seq_iter = self.seqs,
                sample_size = 2,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = ['muscle'])
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertEqual(
                    self.expected[seq1.id][seq2.id],
                    d)
        self.assertEqual(i, 5)

        distance_iter = seqstats.sample_distance_iter(
                seq_iter = self.seqs,
                sample_size = 2,
                per_site = True,
                aligned = False,
                ignore_gaps = True,
                aligner_tools = ['muscle'])
        for i, (seq1, seq2, d, drc) in enumerate(distance_iter):
            self.assertAlmostEqual(
                    self.expected[seq1.id][seq2.id] / 4.0,
                    d)
        self.assertEqual(i, 5)

class DistanceTestCase(unittest.TestCase):
    def test_aligned(self):
        seq1 = SeqRecord(Seq('AC--GTNAC-TYATR'), id='1')
        seq2 = SeqRecord(Seq('ACN-GTAAC--CATT'), id='2')
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = True)
        self.assertEqual(d, 1)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = True)
        self.assertAlmostEqual(dps, 1 / float(15))

        d = seqstats.distance(seq1, seq2, per_site = False, aligned = True,
                ignore_gaps = False)
        self.assertEqual(d, 3)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = True,
                ignore_gaps = False)
        self.assertAlmostEqual(dps, 3 / float(15))

    def test_unaligned(self):
        seq1 = SeqRecord(Seq('AC--GTNAC-TYATR'), id='1')
        seq2 = SeqRecord(Seq('ACN-GTAAC--CATT'), id='2')
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = False,
                ignore_gaps = False,
                aligner_tools = None)
        self.assertEqual(d, 3)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = False,
                ignore_gaps = False,
                aligner_tools = None)
        self.assertAlmostEqual(dps, 3 / float(13))

        seq1 = SeqRecord(Seq('ATCCGT'), id='1')
        seq2 = SeqRecord(Seq('ACCGT'), id='2')
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = False,
                ignore_gaps = True,
                aligner_tools = None)
        self.assertEqual(d, 0)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = False,
                ignore_gaps = True,
                aligner_tools = None)
        self.assertEqual(dps, 0.0)
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = False,
                ignore_gaps = False,
                aligner_tools = None)
        self.assertEqual(d, 1)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = False,
                ignore_gaps = False,
                aligner_tools = None)
        self.assertEqual(dps, 1 / float(6))

    def test_unaligned_mafft(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        seq1 = SeqRecord(Seq('AC--GTNAC-TYATR'), id='1')
        seq2 = SeqRecord(Seq('ACN-GTAAC--CATT'), id='2')
        # 'ACGTNACTYATR'
        # 'ACNGTAACCATT'
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['mafft'])
        self.assertEqual(d, 4)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['mafft'])
        self.assertAlmostEqual(dps, 4 / float(12))

        seq1 = SeqRecord(Seq('ATCCGT'), id='1')
        seq2 = SeqRecord(Seq('ACCGT'), id='2')
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = False,
                ignore_gaps = True,
                aligner_tools = ['mafft'])
        self.assertEqual(d, 0)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = False,
                ignore_gaps = True,
                aligner_tools = ['mafft'])
        self.assertEqual(dps, 0.0)
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['mafft'])
        self.assertEqual(d, 1)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['mafft'])
        self.assertEqual(dps, 1 / float(6))

    def test_unaligned_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        seq1 = SeqRecord(Seq('AC--GTNAC-TYATR'), id='1')
        seq2 = SeqRecord(Seq('ACN-GTAAC--CATT'), id='2')
        # 'ACGTNACTYATR'
        # 'ACNGTAACCATT'
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['muscle'])
        self.assertEqual(d, 4)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['muscle'])
        self.assertAlmostEqual(dps, 4 / float(12))

        seq1 = SeqRecord(Seq('ATCCGT'), id='1')
        seq2 = SeqRecord(Seq('ACCGT'), id='2')
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = False,
                ignore_gaps = True,
                aligner_tools = ['muscle'])
        self.assertEqual(d, 0)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = False,
                ignore_gaps = True,
                aligner_tools = ['muscle'])
        self.assertEqual(dps, 0.0)
        d = seqstats.distance(seq1, seq2, per_site = False, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['muscle'])
        self.assertEqual(d, 1)
        dps = seqstats.distance(seq1, seq2, per_site = True, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['muscle'])
        self.assertEqual(dps, 1 / float(6))

class GetDifferencesTestCase(unittest.TestCase):
    def test_align_error(self):
        seq1 = SeqRecord(Seq('ACGT'), id='1')
        seq2 = SeqRecord(Seq('TACGT'), id='2')
        self.assertRaises(seqstats.AlignmentError, seqstats.get_differences,
                seq1, seq2,
                True)
        diffs, l = seqstats.get_differences(seq1, seq2)
        self.assertEqual(diffs, {})
        self.assertEqual(l, 5)
        expected = {0: ('-', 'T')}
        diffs, l = seqstats.get_differences(seq1, seq2, ignore_gaps = False)
        self.assertEqual(diffs, expected)
        self.assertEqual(l, 5)

    def test_aligned(self):
        seq1 = SeqRecord(Seq('AC--GTNAC-TYATR'), id='1')
        seq2 = SeqRecord(Seq('ACN-GTAAC--CATT'), id='2')
        e = {14:('R', 'T')}
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = True)
        self.assertEqual(diffs, e)
        self.assertEqual(l, 15)
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = True,
                ignore_gaps = False)
        e = {14: ('R', 'T'),
             2:  ('-', 'N'),
             10: ('T', '-')}
        self.assertEqual(diffs, e)
        self.assertEqual(l, 15)

    def test_unaligned(self):
        seq1 = SeqRecord(Seq('AC--GTNAC-TYATR'), id='1')
        seq2 = SeqRecord(Seq('ACN-GTAAC--CATT'), id='2')
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = False,
                ignore_gaps = False,
                aligner_tools = None)
        e = {12: ('R', 'T'),
             2:  ('-', 'N'),
             8: ('T', '-')}
        self.assertEqual(diffs, e)
        self.assertEqual(l, 13)

        seq1 = SeqRecord(Seq('ATCCGT'), id='1')
        seq2 = SeqRecord(Seq('ACCGT'), id='2')
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = False,
                ignore_gaps = True,
                aligner_tools = None)
        self.assertEqual(diffs, {})
        self.assertEqual(l, 6)
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = False,
                ignore_gaps = False,
                aligner_tools = None)
        self.assertEqual(diffs, {1: ('T', '-')})
        self.assertEqual(l, 6)

    def test_unaligned_mafft(self):
        if not functions.which('mafft'):
            _LOG.warning('mafft not found... skipping tests.')
            return
        seq1 = SeqRecord(Seq('AC--GTNAC-TYATR'), id='1')
        seq2 = SeqRecord(Seq('ACN-GTAAC--CATT'), id='2')
        # 'ACGTNACTYATR'
        # 'ACNGTAACCATT'
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['mafft'])
        e = {11: ('R', 'T'),
             3:  ('T', 'G'),
             6:  ('C', 'A'),
             7:  ('T', 'C')}
        self.assertEqual(diffs, e)
        self.assertEqual(l, 12)

        seq1 = SeqRecord(Seq('ATCCGT'), id='1')
        seq2 = SeqRecord(Seq('ACCGT'), id='2')
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = False,
                ignore_gaps = True,
                aligner_tools = ['mafft'])
        self.assertEqual(diffs, {})
        self.assertEqual(l, 6)
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['mafft'])
        self.assertEqual(diffs, {1: ('T', '-')})
        self.assertEqual(l, 6)

    def test_unaligned_muscle(self):
        if not functions.which('muscle'):
            _LOG.warning('muscle not found... skipping tests.')
            return
        seq1 = SeqRecord(Seq('AC--GTNAC-TYATR'), id='1')
        seq2 = SeqRecord(Seq('ACN-GTAAC--CATT'), id='2')
        # 'ACGTNACTYATR'
        # 'ACNGTAACCATT'
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['muscle'])
        e = {11: ('R', 'T'),
             3:  ('T', 'G'),
             6:  ('C', 'A'),
             7:  ('T', 'C')}
        self.assertEqual(diffs, e)
        self.assertEqual(l, 12)

        seq1 = SeqRecord(Seq('ATCCGT'), id='1')
        seq2 = SeqRecord(Seq('ACCGT'), id='2')
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = False,
                ignore_gaps = True,
                aligner_tools = ['muscle'])
        self.assertEqual(diffs, {})
        self.assertEqual(l, 6)
        diffs, l = seqstats.get_differences(seq1, seq2, aligned = False,
                ignore_gaps = False,
                aligner_tools = ['muscle'])
        self.assertEqual(diffs, {1: ('T', '-')})
        self.assertEqual(l, 6)

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
        self.assertRaises(seqstats.AlignmentError, seqstats.column_frequencies,
                self.unaligned)
    
    def test_unambiguous_bases(self):
        for i, char in enumerate(['A', 'C', 'G', 'T']):
            expected_freqs = [0.0] * 5
            expected_freqs[i] = 1.0
            col_freqs, seqs = seqstats.column_frequencies(
                    self.simple_alignment,
                    character_list=[char])
            self.assertEqual(
                    col_freqs,
                    expected_freqs)
            self.assertSameSequences(self.simple_alignment, seqs)

    def test_gap_chars(self):
        expected_freqs = ([0.0] * 4) + [1.0]
        col_freqs, seqs = seqstats.column_frequencies(
                self.simple_alignment)
        self.assertEqual(
                col_freqs,
                expected_freqs)
        self.assertSameSequences(self.simple_alignment, seqs)
        expected_freqs = ([0.0] * 4) + [3/float(5)]
        col_freqs, seqs = seqstats.column_frequencies(
                self.simple_alignment,
                character_list = ['?'])
        self.assertAlmostEqual(
                col_freqs,
                expected_freqs)
        self.assertSameSequences(self.simple_alignment, seqs)
    
    def test_missing_character(self):
        expected_freqs = [0.0] * 5
        col_freqs, seqs = seqstats.column_frequencies(
                self.simple_alignment,
                character_list = ['x'])
        self.assertAlmostEqual(
                col_freqs,
                expected_freqs)
        self.assertSameSequences(self.simple_alignment, seqs)

if __name__ == '__main__':
    unittest.main()

