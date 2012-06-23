#! /usr/bin/env python

import os
import sys
import unittest
import itertools

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

from seqsift.digest import *
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class FragmentTestCase(unittest.TestCase):
    def test_init_str(self):
        f = Fragment('ACGT',
                start_site=10,
                end_site=13,
                overhang=2,
                five_prime_terminus=True,
                three_prime_terminus=False,
                id='test')
        self.assertIsInstance(f, Fragment)
        self.assertEqual(f.start, 10)
        self.assertEqual(f.end, 13)
        self.assertEqual(f.overhang, 2)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(f.length, 4)
        self.assertEqual(str(f.seq), 'ACGT')
        self.assertEqual(f.id, 'test')
        self.assertEqual(f.name, '<unknown name>')

    def test_init_seq(self):
        f = Fragment(Seq('ACGT'),
                start_site=10,
                end_site=13,
                overhang=2,
                five_prime_terminus=True,
                three_prime_terminus=False,
                id='test')
        self.assertIsInstance(f, Fragment)
        self.assertEqual(f.start, 10)
        self.assertEqual(f.end, 13)
        self.assertEqual(f.overhang, 2)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(f.length, 4)
        self.assertEqual(str(f.seq), 'ACGT')
        self.assertEqual(f.id, 'test')
        self.assertEqual(f.name, '<unknown name>')

    def test_init_seq_record(self):
        f = Fragment(SeqRecord(Seq('ACGT'), id='test'),
                start_site=10,
                end_site=13,
                overhang=2,
                five_prime_terminus=True,
                three_prime_terminus=False)
        self.assertIsInstance(f, Fragment)
        self.assertEqual(f.start, 10)
        self.assertEqual(f.end, 13)
        self.assertEqual(f.overhang, 2)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(f.length, 4)
        self.assertEqual(str(f.seq), 'ACGT')
        self.assertEqual(f.id, 'test')
        self.assertEqual(f.name, '<unknown name>')

    def test_init_site_index_error(self):
        self.assertRaises(InvalidFragmentError, Fragment, 'ACGT', 2, 1, 0)

    def test_init_length_error(self):
        self.assertRaises(InvalidFragmentError, Fragment, 'ACGT', 1, 10, 0)

class RecognitionSeqTestCase(unittest.TestCase):
    def setUp(self):
        self.ecori = SeqRecord(Seq('GAATTC', IUPAC.unambiguous_dna),
            id = 'ecori',
            name = 'EcoRI',
            description = 'Recognition sequence of EcoRI restriction enzyme')
        self.sphi = SeqRecord(Seq('GCATGC', IUPAC.unambiguous_dna),
            id = 'sphi',
            name = 'SphI',
            description = 'Recognition sequence of SphI restriction enzyme')

    def test_init_str(self):
        rs = RecognitionSeq(str(self.ecori.seq), 1)
        self.assertIsInstance(rs, RecognitionSeq)
        self.assertIsInstance(rs.seq, Seq)

    def test_init_seq(self):
        rs = RecognitionSeq(self.ecori.seq, 1)
        self.assertIsInstance(rs, RecognitionSeq)
        self.assertIsInstance(rs.seq, Seq)

    def test_init_seq_record(self):
        rs = RecognitionSeq(self.ecori, 1)
        self.assertIsInstance(rs, RecognitionSeq)
        self.assertIsInstance(rs.seq, Seq)

    def test_init_invalid_base(self):
        self.assertRaises(InvalidRecognitionSeqError,
                RecognitionSeq, 'AGCU', 0)

    def test_init_invalid_cut_site(self):
        self.assertRaises(InvalidRecognitionSeqError,
                RecognitionSeq, 'AGCU', -1)
        self.assertRaises(InvalidRecognitionSeqError,
                RecognitionSeq, 'AGCU', 5)

    def test_overhang(self):
        rs = RecognitionSeq(self.ecori, 0)
        self.assertEqual(rs.overhang, 6)
        rs = RecognitionSeq(self.ecori, 6)
        self.assertEqual(rs.overhang, -6)
        rs = RecognitionSeq(self.ecori, 3)
        self.assertEqual(rs.overhang, 0)

    def test_palindrome(self):
        rs = RecognitionSeq(self.ecori, 1)
        self.assertTrue(rs.palindrome)
        rs = RecognitionSeq('GATATC', 3)
        self.assertTrue(rs.palindrome)
        rs = RecognitionSeq('TCGA', 1)
        self.assertTrue(rs.palindrome)
        rs = RecognitionSeq('TAAT', 2)
        self.assertFalse(rs.palindrome)

    def test_potential_overlap(self):
        for i in range(len(self.ecori)+1):
            rs = RecognitionSeq(self.ecori, i)
            self.assertFalse(rs.potential_overlap)
        for i in range(len(self.sphi)+1):
            rs = RecognitionSeq(self.sphi, i)
            if i < 2 or i > 4:
                self.assertTrue(rs.potential_overlap)
            else:
                self.assertFalse(rs.potential_overlap)

    def test_short(self):
        rs = RecognitionSeq(self.ecori, 1)
        mol = self.ecori[2:]
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 1)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 4)
        self.assertEqual(f.length, 4)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq[2:]))

    def test_no_cut(self):
        rs = RecognitionSeq(self.ecori, 1)
        mol = SeqRecord(Seq('AAAAAAAAAA'), id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 1)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 10)
        self.assertEqual(f.length, 10)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), 'AAAAAAAAAA')

    def test_five_prime_min_case(self):
        rs = RecognitionSeq(self.ecori, 0)
        fragments = list(rs.digest(self.ecori))
        self.assertEqual(len(fragments), 1)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 6)
        self.assertEqual(f.length, 6)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq))
    
    def test_three_prime_min_case(self):
        rs = RecognitionSeq(self.ecori, 6)
        fragments = list(rs.digest(self.ecori))
        self.assertEqual(len(fragments), 1)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 6)
        self.assertEqual(f.length, 6)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq))

    def test_blunt_min_case(self):
        rs = RecognitionSeq(self.ecori, 3)
        fragments = list(rs.digest(self.ecori))
        self.assertEqual(len(fragments), 2)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 3)
        self.assertEqual(f.length, 3)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq)[:3])
        f = fragments[1]
        self.assertEqual(f.start, 4)
        self.assertEqual(f.end, 6)
        self.assertEqual(f.length, 3)
        self.assertEqual(f.overhang, 0)
        self.assertFalse(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq)[3:])

    def test_five_prime_head_cut_case(self):
        rs = RecognitionSeq(self.ecori, 0)
        mol = SeqRecord(Seq(str(self.ecori.seq) + 'AAAAAAAAAA'), id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 1)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 16)
        self.assertEqual(f.length, 16)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(mol.seq))

    def test_three_prime_head_cut_case(self):
        rs = RecognitionSeq(self.ecori, 6)
        mol = SeqRecord(Seq(str(self.ecori.seq) + 'AAAAAAAAAA'), id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 2)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 6)
        self.assertEqual(f.length, 6)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq))
        f = fragments[1]
        self.assertEqual(f.start, 7)
        self.assertEqual(f.end, 16)
        self.assertEqual(f.length, 10)
        self.assertEqual(f.overhang, 6)
        self.assertFalse(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), 'AAAAAAAAAA')

    def test_five_prime_tail_cut_case(self):
        rs = RecognitionSeq(self.ecori, 0)
        mol = SeqRecord(Seq( 'AAAAAAAAAA' + str(self.ecori.seq)), id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 2)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 10)
        self.assertEqual(f.length, 10)
        self.assertEqual(f.overhang, 6)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(str(f.seq), 'AAAAAAAAAA')
        f = fragments[1]
        self.assertEqual(f.start, 11)
        self.assertEqual(f.end, 16)
        self.assertEqual(f.length, 6)
        self.assertEqual(f.overhang, 0)
        self.assertFalse(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq))

    def test_three_prime_tail_cut_case(self):
        rs = RecognitionSeq(self.ecori, 6)
        mol = SeqRecord(Seq( 'AAAAAAAAAA' + str(self.ecori.seq)), id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 1)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 16)
        self.assertEqual(f.length, 16)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(mol.seq))

    def test_five_prime_repeat_case(self):
        rs = RecognitionSeq(self.ecori, 0)
        mol = SeqRecord(Seq(str(self.ecori.seq) * 3), id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 3)
        for i, f in enumerate(fragments):
            self.assertEqual(f.start, i*(len(rs))+1)
            self.assertEqual(f.end, i*(len(rs))+6)
            self.assertEqual(f.length, 6)
            if i == 2:
                self.assertEqual(f.overhang, 0)
            else:
                self.assertEqual(f.overhang, 6)
            if i == 0:
                self.assertTrue(f.five_prime_terminus)
            else:
                self.assertFalse(f.five_prime_terminus)
            if i == 2:
                self.assertTrue(f.three_prime_terminus)
            else:
                self.assertFalse(f.three_prime_terminus)
            self.assertEqual(str(f.seq), str(self.ecori.seq))

    def test_three_prime_repeat_case(self):
        rs = RecognitionSeq(self.ecori, 6)
        mol = SeqRecord(Seq(str(self.ecori.seq) * 3), id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 3)
        for i, f in enumerate(fragments):
            self.assertEqual(f.start, i*(len(rs))+1)
            self.assertEqual(f.end, i*(len(rs))+6)
            self.assertEqual(f.length, 6)
            if i == 0:
                self.assertEqual(f.overhang, 0)
            else:
                self.assertEqual(f.overhang, 6)
            if i == 0:
                self.assertTrue(f.five_prime_terminus)
            else:
                self.assertFalse(f.five_prime_terminus)
            if i == 2:
                self.assertTrue(f.three_prime_terminus)
            else:
                self.assertFalse(f.three_prime_terminus)
            self.assertEqual(str(f.seq), str(self.ecori.seq))

    def test_blunt_repeat_case(self):
        rs = RecognitionSeq(self.ecori, 3)
        mol = SeqRecord(Seq(str(self.ecori.seq) * 3), id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 4)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 3)
        self.assertEqual(f.length, 3)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq)[:3])
        for i in range(1,3):
            f = fragments[i]
            self.assertEqual(f.start, i*len(self.ecori)-2)
            self.assertEqual(f.end, i*len(self.ecori)+3)
            self.assertEqual(f.length, 6)
            self.assertEqual(f.overhang, 0)
            self.assertFalse(f.five_prime_terminus)
            self.assertFalse(f.three_prime_terminus)
            self.assertEqual(str(f.seq), str(self.ecori.seq)[3:] + 
                    str(self.ecori.seq)[:3])
        f = fragments[3]
        self.assertEqual(f.start, 16)
        self.assertEqual(f.end, 18)
        self.assertEqual(f.length, 3)
        self.assertEqual(f.overhang, 0)
        self.assertFalse(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq)[3:])

    def test_five_prime_one_cut_case(self):
        rs = RecognitionSeq(self.ecori, 0)
        mol = SeqRecord(
                Seq( 'AAAAAAAAAA' + str(self.ecori.seq) + 'TTTTTTTTTT'),
                id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 2)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 10)
        self.assertEqual(f.length, 10)
        self.assertEqual(f.overhang, 6)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(str(f.seq), 'AAAAAAAAAA')
        f = fragments[1]
        self.assertEqual(f.start, 11)
        self.assertEqual(f.end, 26)
        self.assertEqual(f.length, 16)
        self.assertEqual(f.overhang, 0)
        self.assertFalse(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq) + 'TTTTTTTTTT')

    def test_three_prime_one_cut_case(self):
        rs = RecognitionSeq(self.ecori, 6)
        mol = SeqRecord(
                Seq( 'AAAAAAAAAA' + str(self.ecori.seq) + 'TTTTTTTTTT'),
                id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 2)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 16)
        self.assertEqual(f.length, 16)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(str(f.seq), 'AAAAAAAAAA' + str(self.ecori.seq))
        f = fragments[1]
        self.assertEqual(f.start, 17)
        self.assertEqual(f.end, 26)
        self.assertEqual(f.length, 10)
        self.assertEqual(f.overhang, 6)
        self.assertFalse(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), 'TTTTTTTTTT')

    def test_blunt_one_cut_case(self):
        rs = RecognitionSeq(self.ecori, 3)
        mol = SeqRecord(
                Seq( 'AAAAAAAAAA' + str(self.ecori.seq) + 'TTTTTTTTTT'),
                id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 2)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 13)
        self.assertEqual(f.length, 13)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(str(f.seq), 'AAAAAAAAAA' + str(self.ecori.seq)[:3])
        f = fragments[1]
        self.assertEqual(f.start, 14)
        self.assertEqual(f.end, 26)
        self.assertEqual(f.length, 13)
        self.assertEqual(f.overhang, 0)
        self.assertFalse(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq)[3:] + 'TTTTTTTTTT')

    def test_cut_site_two_one_cut_case(self):
        rs = RecognitionSeq(self.ecori, 2)
        mol = SeqRecord(
                Seq( 'AAAAAAAAAA' + str(self.ecori.seq) + 'TTTTTTTTTT'),
                id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 2)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 12)
        self.assertEqual(f.length, 12)
        self.assertEqual(f.overhang, 2)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(str(f.seq), 'AAAAAAAAAA' + str(self.ecori.seq)[:2])
        f = fragments[1]
        self.assertEqual(f.start, 13)
        self.assertEqual(f.end, 26)
        self.assertEqual(f.length, 14)
        self.assertEqual(f.overhang, 0)
        self.assertFalse(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq)[2:] + 'TTTTTTTTTT')

    def test_cut_site_four_one_cut_case(self):
        rs = RecognitionSeq(self.ecori, 4)
        mol = SeqRecord(
                Seq( 'AAAAAAAAAA' + str(self.ecori.seq) + 'TTTTTTTTTT'),
                id='test')
        fragments = list(rs.digest(mol))
        self.assertEqual(len(fragments), 2)
        f = fragments[0]
        self.assertEqual(f.start, 1)
        self.assertEqual(f.end, 14)
        self.assertEqual(f.length, 14)
        self.assertEqual(f.overhang, 0)
        self.assertTrue(f.five_prime_terminus)
        self.assertFalse(f.three_prime_terminus)
        self.assertEqual(str(f.seq), 'AAAAAAAAAA' + str(self.ecori.seq)[:4])
        f = fragments[1]
        self.assertEqual(f.start, 15)
        self.assertEqual(f.end, 26)
        self.assertEqual(f.length, 12)
        self.assertEqual(f.overhang, 2)
        self.assertFalse(f.five_prime_terminus)
        self.assertTrue(f.three_prime_terminus)
        self.assertEqual(str(f.seq), str(self.ecori.seq)[4:] + 'TTTTTTTTTT')

    def test_simble_gb_seq(self):
        rs = RecognitionSeq('TAG', 3)
        fp = package_paths.data_path('JF314863-JF314866.gb')
        seqs = SeqIO.parse(fp, format='gb', alphabet=IUPAC.ambiguous_dna)
        s = seqs.next()
        self.assertEqual(s.name, 'JF314863')
        fragments = list(rs.digest(s))
        self.assertEqual(len(fragments), 6)
        lengths = [117, 172, 62, 10, 10, 102]
        for i in range(len(fragments)):
            f = fragments[i]
            self.assertIsInstance(f, Fragment)
            self.assertEqual(len(f), lengths[i])
            

class DigestSummaryTestCase(unittest.TestCase):
    def setUp(self):
        self.ecori = RecognitionSeq(Seq('GAATTC', IUPAC.unambiguous_dna),
            cut_site = 3,
            id = 'ecori',
            name = 'EcoRI',
            description = 'Recognition sequence of EcoRI restriction enzyme')
    
    def test_basic(self):
        mol = SeqRecord(Seq(str(self.ecori.seq) * 3), id='test')
        ds = DigestSummary(self.ecori, mol)
        self.assertIsInstance(ds, DigestSummary)
        self.assertEqual(ds.recognition_seq, str(self.ecori.seq))
        self.assertEqual(ds.molecule_id, mol.id)
        self.assertEqual(ds.molecule_name, mol.name)
        self.assertEqual(ds.molecule_description, mol.description)
        self.assertIsInstance(ds.length_distribution, dict)
        # self.assertEqual(ds.length_distribution, {3: 2, 6: 2})
        self.assertEqual(ds.length_distribution, {6: 2})
        self.assertEqual(ds.molecule_length, len(mol))

    def test_simble_gb_seq(self):
        rs = RecognitionSeq('TAG', 3)
        fp = package_paths.data_path('JF314863-JF314866.gb')
        seqs = SeqIO.parse(fp, format='gb', alphabet=IUPAC.ambiguous_dna)
        s = seqs.next()
        self.assertEqual(s.name, 'JF314863')
        ds = DigestSummary(rs, s)
        self.assertIsInstance(ds, DigestSummary)
        self.assertEqual(ds.recognition_seq, str(rs.seq))
        self.assertEqual(ds.molecule_id, s.id)
        self.assertEqual(ds.molecule_name, s.name)
        self.assertEqual(ds.molecule_description, s.description)
        self.assertIsInstance(ds.length_distribution, dict)
        self.assertEqual(ds.length_distribution, {
                13: 2,
                65: 1,
                175: 1,})
        self.assertEqual(ds.molecule_length, len(s))

if __name__ == '__main__':
    unittest.main()
