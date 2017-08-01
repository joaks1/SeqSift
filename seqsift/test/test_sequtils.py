#! /usr/bin/env python

import os
import sys
import types
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

import seqsift
from seqsift.seqops import sequtils
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class SeqBatchIterTestCase(SeqSiftTestCase):
    def test_even(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTGA'), id='1'),
                SeqRecord(Seq('ATGACCAACTGT'), id='2'),
                SeqRecord(Seq('ATGACCAACTGC'), id='3'),
                SeqRecord(Seq('ATGACCAACTGG'), id='4'),
                SeqRecord(Seq('ATGACCAACTAT'), id='5'),
                SeqRecord(Seq('ATGACCAACTAC'), id='6'),
                ]
        batch_iter = sequtils.seq_batch_iter(seqs, 2)
        self.assertIsInstance(batch_iter, types.GeneratorType)
        batch_iter = list(batch_iter)
        self.assertEqual(len(batch_iter), 3)
        new_seqs = []
        for seq_iter in batch_iter:
            self.assertIsInstance(seq_iter, seqsift.utils.dataio.BufferedIter)
            for s in seq_iter:
                new_seqs.append(s)
        self.assertSameData(seqs, new_seqs)

    def test_odd(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTGA'), id='1'),
                SeqRecord(Seq('ATGACCAACTGT'), id='2'),
                SeqRecord(Seq('ATGACCAACTGC'), id='3'),
                SeqRecord(Seq('ATGACCAACTGG'), id='4'),
                SeqRecord(Seq('ATGACCAACTAT'), id='5'),
                SeqRecord(Seq('ATGACCAACTAC'), id='6'),
                SeqRecord(Seq('TTGACCAACTAC'), id='7'),
                ]
        batch_iter = sequtils.seq_batch_iter(seqs, 3)
        self.assertIsInstance(batch_iter, types.GeneratorType)
        batch_iter = list(batch_iter)
        self.assertEqual(len(batch_iter), 3)
        new_seqs = []
        for seq_iter in batch_iter:
            self.assertIsInstance(seq_iter, seqsift.utils.dataio.BufferedIter)
            for s in seq_iter:
                new_seqs.append(s)
        self.assertSameData(seqs, new_seqs)

        batch_iter = sequtils.seq_batch_iter(seqs, 4)
        self.assertIsInstance(batch_iter, types.GeneratorType)
        batch_iter = list(batch_iter)
        self.assertEqual(len(batch_iter), 2)
        new_seqs = []
        for seq_iter in batch_iter:
            self.assertIsInstance(seq_iter, seqsift.utils.dataio.BufferedIter)
            for s in seq_iter:
                new_seqs.append(s)
        self.assertSameData(seqs, new_seqs)

    def test_one_even(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTGA'), id='1'),
                SeqRecord(Seq('ATGACCAACTGT'), id='2'),
                SeqRecord(Seq('ATGACCAACTGC'), id='3'),
                SeqRecord(Seq('ATGACCAACTGG'), id='4'),
                SeqRecord(Seq('ATGACCAACTAT'), id='5'),
                SeqRecord(Seq('ATGACCAACTAC'), id='6'),
                ]
        batch_iter = sequtils.seq_batch_iter(seqs, 6)
        self.assertIsInstance(batch_iter, types.GeneratorType)
        batch_iter = list(batch_iter)
        self.assertEqual(len(batch_iter), 1)
        new_seqs = []
        for seq_iter in batch_iter:
            self.assertIsInstance(seq_iter, seqsift.utils.dataio.BufferedIter)
            for s in seq_iter:
                new_seqs.append(s)
        self.assertSameData(seqs, new_seqs)

    def test_one_odd(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTGA'), id='1'),
                SeqRecord(Seq('ATGACCAACTGT'), id='2'),
                SeqRecord(Seq('ATGACCAACTGC'), id='3'),
                SeqRecord(Seq('ATGACCAACTGG'), id='4'),
                SeqRecord(Seq('ATGACCAACTAT'), id='5'),
                SeqRecord(Seq('ATGACCAACTAC'), id='6'),
                ]
        batch_iter = sequtils.seq_batch_iter(seqs, 10)
        self.assertIsInstance(batch_iter, types.GeneratorType)
        batch_iter = list(batch_iter)
        self.assertEqual(len(batch_iter), 1)
        new_seqs = []
        for seq_iter in batch_iter:
            self.assertIsInstance(seq_iter, seqsift.utils.dataio.BufferedIter)
            for s in seq_iter:
                new_seqs.append(s)
        self.assertSameData(seqs, new_seqs)

class SequencesAreEqualTestCase(SeqSiftTestCase):
    def setUp(self):
        self.seq = SeqIO.read(
                package_paths.data_path('JF314862.gb'),
                format='gb',
                alphabet=IUPAC.ambiguous_dna)

    def test_copy(self):
        seq2 = SeqIO.read(
                package_paths.data_path('JF314862.gb'),
                format='gb',
                alphabet=IUPAC.ambiguous_dna)
        self.assertTrue(sequtils.sequences_are_equal(self.seq, seq2))
        seq2.name += 'a'
        self.assertFalse(sequtils.sequences_are_equal(self.seq, seq2))

class GetLongestReadingFrameTestCase(SeqSiftTestCase):
    def test_cds(self):
        seq = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        exp = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = False,
                require_start_after_stop = False)
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertFalse(lrf is seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

        seq = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        exp = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = False,
                require_start_after_stop = True)
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertFalse(lrf is seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

        seq = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        exp = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = True,
                require_start_after_stop = False)
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertNotEqual(lrf, seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

        seq = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        exp = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = True,
                require_start_after_stop = True)
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertNotEqual(lrf, seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

    def test_no_partial(self):
        seq = SeqRecord(Seq('AAGACCAACTGAATA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = False,
                require_start_after_stop = True)
        self.assertEqual(lrf, [])
        seq = SeqRecord(Seq('ATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1')
        exp = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = False,
                require_start_after_stop = True)
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertFalse(lrf is seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

    def test_partial(self):
        seq = SeqRecord(Seq('AAGACCAACTGAATA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = True,
                require_start_after_stop = True)
        exp = SeqRecord(Seq('AGACCAACTGAATA', IUPAC.ambiguous_dna), id='1')
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertFalse(lrf is seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

    def test_require_start_after_stop(self):
        seq = SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1')
        exp = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = True,
                require_start_after_stop = True)
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertFalse(lrf is seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

        seq = SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1')
        exp = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = False,
                require_start_after_stop = True)
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertNotEqual(lrf, seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

        seq = SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1')
        exp = SeqRecord(Seq('ATAGAAATTGGCCATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = True,
                require_start_after_stop = False)
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertNotEqual(lrf, seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

        seq = SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1')
        exp = SeqRecord(Seq('ATGACCAACTGA', IUPAC.ambiguous_dna), id='1')
        lrf = sequtils.get_longest_reading_frames(seq,
                table = 1,
                allow_partial = False,
                require_start_after_stop = False)
        self.assertEqual(len(lrf), 1)
        lrf = lrf[0]
        self.assertNotEqual(lrf, seq)
        self.assertSameMetadata(lrf, seq)
        self.assertEqual(str(lrf.seq), str(exp.seq))

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
        self.assertFalse(s is self.seq)

    def test_string(self):
        s = sequtils.copy_seq_metadata(self.seq, 'AGCT')
        self.assertSameMetadata(self.seq, s)
        self.assertEqual(str(s.seq), 'AGCT')
        self.assertFalse(s is self.seq)

    def test_seq(self):
        s = sequtils.copy_seq_metadata(self.seq, Seq('AGCT'))
        self.assertSameMetadata(self.seq, s)
        self.assertEqual(str(s.seq), 'AGCT')
        self.assertFalse(s is self.seq)

    def test_seq_record(self):
        s = sequtils.copy_seq_metadata(self.seq, SeqRecord(Seq('AGCT'), id='1'))
        self.assertSameMetadata(self.seq, s)
        self.assertEqual(str(s.seq), 'AGCT')
        self.assertFalse(s is self.seq)

class GetWithoutGapsTestCase(SeqSiftTestCase):
    def test_without_gaps(self):
        s1 = SeqRecord(Seq('CCATG'), id='1')
        s2 = sequtils.get_without_gaps(s1)
        self.assertSameMetadata(s1, s2)
        self.assertEqual(str(s1.seq), str(s2.seq))
        self.assertFalse(s1 is s2)

    def test_without_gaps(self):
        s1 = SeqRecord(Seq('--CC-AT--G---'), id='1')
        s2 = sequtils.get_without_gaps(s1)
        self.assertSameMetadata(s1, s2)
        self.assertEqual(str(s2.seq), 'CCATG')
        self.assertFalse(s1 is s2)

class GetTranslation(SeqSiftTestCase):
    def test_unambiguous_dna(self):
        s1 = SeqRecord(Seq('ATGACCAACTGA', IUPAC.unambiguous_dna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertFalse(s1 is s2)
        self.assertEqual(str(s2.seq), 'MTN*')

    def test_unambiguous_rna(self):
        s1 = SeqRecord(Seq('AUGACCAACUGA', IUPAC.unambiguous_rna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertFalse(s1 is s2)
        self.assertEqual(str(s2.seq), 'MTN*')

    def test_ambiguous_dna(self):
        s1 = SeqRecord(Seq('ATGATRAACTGA', IUPAC.ambiguous_dna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertFalse(s1 is s2)
        self.assertEqual(str(s2.seq), 'MXN*')

    def test_ambiguous_rna(self):
        s1 = SeqRecord(Seq('AUGAURAACUGA', IUPAC.ambiguous_rna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertFalse(s1 is s2)
        self.assertEqual(str(s2.seq), 'MXN*')

    def test_to_stop(self):
        s1 = SeqRecord(Seq('ATGACCAACTGA', IUPAC.unambiguous_dna), id='1')
        s2 = sequtils.get_translation(s1, to_stop = True)
        self.assertSameMetadata(s1, s2)
        self.assertFalse(s1 is s2)
        self.assertEqual(str(s2.seq), 'MTN')

    def test_extra_base(self):
        s1 = SeqRecord(Seq('ATGACCAACTGAATAA', IUPAC.unambiguous_dna), id='1')
        s2 = sequtils.get_translation(s1)
        self.assertSameMetadata(s1, s2)
        self.assertFalse(s1 is s2)
        self.assertEqual(str(s2.seq), 'MTN*I')

if __name__ == '__main__':
    unittest.main()
