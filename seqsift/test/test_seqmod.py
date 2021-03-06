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

class LongestReadingFramesTestCase(SeqSiftTestCase):
    def test_simple(self):
        seqs = [
                SeqRecord(Seq('--TAGAT-AGAT--AGA-ATTGG-CCA-TGAC-CAAC-TGA-ATA'), id='1'),
                SeqRecord(Seq('T-AGA--TAGATAG-AAATTG-GCCAT--GAC-C-AA-GTGA-ATA-'), id='2'),
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACAAACTGAATA'), id='3'),
                ]
        exp = [
                SeqRecord(Seq('ATGACCAACTGA'), id='1'),
                SeqRecord(Seq('ATGACCAAGTGA'), id='2'),
                SeqRecord(Seq('ATGACAAACTGA'), id='3'),
                ]
        new_seqs = seqmod.longest_reading_frames(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

class TranslateLongestReadingFramesTestCase(SeqSiftTestCase):
    def test_simple(self):
        seqs = [
                SeqRecord(Seq('--TAGAT-AGAT--AGA-ATTGG-CCA-TGAC-CAAC-TGA-ATA'), id='1'),
                SeqRecord(Seq('T-AGA--TAGATAG-AAATTG-GCCAT--GAC-C-AA-GTGA-ATA-'), id='2'),
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACAAACTGAATA'), id='3'),
                ]
        exp = [
                SeqRecord(Seq('MTN*'), id='1'),
                SeqRecord(Seq('MTK*'), id='2'),
                SeqRecord(Seq('MTN*'), id='3'),
                ]
        new_seqs = seqmod.translate_longest_reading_frames(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

class LongestReadingFramesTestCase(SeqSiftTestCase):
    def test_simple(self):
        seqs = [
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACCAACTGAATA'), id='1'),
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACCAAGTGAATA'), id='2'),
                SeqRecord(Seq('TAGATAGATAGAAATTGGCCATGACAAACTGAATA'), id='3'),
                ]
        exp = [
                SeqRecord(Seq('ATGACCAACTGA'), id='1'),
                SeqRecord(Seq('ATGACCAAGTGA'), id='2'),
                SeqRecord(Seq('ATGACAAACTGA'), id='3'),
                ]
        new_seqs = seqmod.longest_reading_frames(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

class ReverseComplementTestCase(SeqSiftTestCase):
    def test_simple(self):
        seqs = [
                SeqRecord(Seq('AAATTACGGG'), id='1'),
                SeqRecord(Seq('TAATTACGGG'), id='2'),
                SeqRecord(Seq('CAATTACGGG'), id='3'),
                SeqRecord(Seq('GAATTACGGG'), id='4'),
                ]
        exp = [
                SeqRecord(Seq('CCCGTAATTT'), id='1'),
                SeqRecord(Seq('CCCGTAATTA'), id='2'),
                SeqRecord(Seq('CCCGTAATTG'), id='3'),
                SeqRecord(Seq('CCCGTAATTC'), id='4'),
                ]
        new_seqs = seqmod.reverse_complement(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

class ReverseComplementToFirstSeqTestCase(SeqSiftTestCase):
    def test_simple(self):
        seqs = [
                SeqRecord(Seq('ATTACGT'), id='1'),
                SeqRecord(Seq('ATTACGT'), id='2'),
                SeqRecord(Seq('ACGTAAT'), id='3'),
                SeqRecord(Seq('ATTACGT'), id='4'),
                ]
        exp = [
                SeqRecord(Seq('ATTACGT'), id='1'),
                SeqRecord(Seq('ATTACGT'), id='2'),
                SeqRecord(Seq('ATTACGT'), id='3'),
                SeqRecord(Seq('ATTACGT'), id='4'),
                ]
        new_seqs = seqmod.reverse_complement_to_first_seq(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)
    
    def test_ambiguity(self):
        seqs = [
                SeqRecord(Seq('AAGG?CCGA', alphabet=IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('AAGG?CCGC', alphabet=IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('AAGG?CCGG', alphabet=IUPAC.ambiguous_dna), id='3'),
                SeqRecord(Seq('ACGG?CCTT', alphabet=IUPAC.ambiguous_dna), id='4'),
                SeqRecord(Seq('AAGG?CCGA', alphabet=IUPAC.ambiguous_dna), id='5'),
                ]
        exp = [
                SeqRecord(Seq('AAGG?CCGA', alphabet=IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('AAGG?CCGC', alphabet=IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('AAGG?CCGG', alphabet=IUPAC.ambiguous_dna), id='3'),
                SeqRecord(Seq('AAGG?CCGT', alphabet=IUPAC.ambiguous_dna), id='4'),
                SeqRecord(Seq('AAGG?CCGA', alphabet=IUPAC.ambiguous_dna), id='5'),
                ]
        new_seqs = seqmod.reverse_complement_to_first_seq(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

class ReverseComplementToLongestReadingFrame(SeqSiftTestCase):
    def test_simple(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTCACTA', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ATGACCAACTCACAC', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('TAGTAAGTTGGTCAT', IUPAC.ambiguous_dna), id='3'),
                ]
        exp = [
                SeqRecord(Seq('ATGACCAACTCACTA', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ATGACCAACTCACAC', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ATGACCAACTTACTA', IUPAC.ambiguous_dna), id='3'),
                ]
        new_seqs = seqmod.reverse_complement_to_longest_reading_frame(seqs)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

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

class DiceTestCase(SeqSiftTestCase):
    def test_one_region(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ATGACCAACTGACAC', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ATGACCAACTGACCC', IUPAC.ambiguous_dna), id='3'),
                ]
        exp = [
                SeqRecord(Seq('CAACTGAA', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('CAACTGAC', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('CAACTGAC', IUPAC.ambiguous_dna), id='3'),
                ]
        region = [(5, 13)]
        new_seqs = seqmod.dice(seqs, region)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

        region = [(0, len(seqs[0]))]
        new_seqs = seqmod.dice(seqs, region)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, seqs)

        region = [(5, -2)]
        new_seqs = seqmod.dice(seqs, region)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

    def test_two_regions(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ATGACCAACTGACAC', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ATGACCAACTGACCC', IUPAC.ambiguous_dna), id='3'),
                ]
        exp = [
                SeqRecord(Seq('GACCGAAT', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('GACCGACA', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('GACCGACC', IUPAC.ambiguous_dna), id='3'),
                ]
        region = [(2, 6), (10, 14)]
        new_seqs = seqmod.dice(seqs, region)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

    def test_two_regions_overlap(self):
        seqs = [
                SeqRecord(Seq('ATGACCAACTGAATA', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ATGACCAACTGACAC', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ATGACCAACTGACCC', IUPAC.ambiguous_dna), id='3'),
                ]
        exp = [
                SeqRecord(Seq('CCAAAACTGA', IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('CCAAAACTGA', IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('CCAAAACTGA', IUPAC.ambiguous_dna), id='3'),
                ]
        region = [(4, 8), (6, 12)]
        new_seqs = seqmod.dice(seqs, region)
        self.assertIsInstance(new_seqs, types.GeneratorType)
        self.assertSameData(new_seqs, exp)

if __name__ == '__main__':
    unittest.main()

