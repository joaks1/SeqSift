#! /usr/bin/env python

import os
import sys
import unittest
import itertools

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from seqsift.utils.dataio import *
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class GetStateAlphabetTestCase(unittest.TestCase):
    def test_error_handling(self):
        self.assertRaises(ValueError, get_state_alphabet, 'bogus')
    def test_ambiguous_dna(self):
        alphabet = get_state_alphabet('dna', ambiguities=True)
        self.assertIsInstance(alphabet, IUPAC.IUPACAmbiguousDNA)
    def test_unambiguous_dna(self):
        alphabet = get_state_alphabet('dna', ambiguities=False)
        self.assertIsInstance(alphabet, IUPAC.IUPACUnambiguousDNA)
    def test_ambiguous_rna(self):
        alphabet = get_state_alphabet('rna', ambiguities=True)
        self.assertIsInstance(alphabet, IUPAC.IUPACAmbiguousRNA)
    def test_unambiguous_rna(self):
        alphabet = get_state_alphabet('rna', ambiguities=False)
        self.assertIsInstance(alphabet, IUPAC.IUPACUnambiguousRNA)
    def test_ambiguous_protein(self):
        for arg in ['protein', 'aa']:
            alphabet = get_state_alphabet(arg, ambiguities=True)
            self.assertIsInstance(alphabet, IUPAC.ExtendedIUPACProtein)
    def test_unambiguous_protein(self):
        for arg in ['protein', 'aa']:
            alphabet = get_state_alphabet(arg, ambiguities=False)
            self.assertIsInstance(alphabet, IUPAC.IUPACProtein)

class BufferedIterTestCase(SeqSiftTestCase):
    def test_simple(self):
        objs = [SeqRecord(Seq('AAAA'), id='a'),
                SeqRecord(Seq('CCCC'), id='c'),
                SeqRecord(Seq('GGGG'), id='g'),
                SeqRecord(Seq('TTTT'), id='t')]
        buffered_iter = BufferedIter(objs)
        for i in range(3):
            self.assertSameIDs(objs, buffered_iter)
            self.assertSameSequences(objs, buffered_iter)

class GetTmpHandleTestCase(unittest.TestCase):
    def setUp(self):
        self.path = package_paths.data_path('primates.partitions.txt')
        self.line = 'another line\n'

    def test_empty_handle(self):
        tmp = get_tmp_handle()
        self.assertTrue(hasattr(tmp, 'close'))
        self.assertEqual(tmp.read(), '')
        tmp.close()

    def test_handle_bottom(self):
        tmp = get_tmp_handle(self.path, rewind=False)
        self.assertTrue(hasattr(tmp, 'close'))
        self.assertEqual(tmp.read(), '')
        tmp.write(self.line)
        tmp.seek(0)
        self.assertEqual(
                tmp.read(),
                open(self.path, 'rU').read() + self.line)
        tmp.close()

    def test_handle_top(self):
        tmp = get_tmp_handle(self.path, rewind=True)
        self.assertTrue(hasattr(tmp, 'close'))
        self.assertEqual(tmp.read(), open(self.path, 'rU').read())

class ReadSeqTestCase(SeqSiftTestCase):
    def setUp(self):
        self.mkTestDir()
        self.gb_path = package_paths.data_path('JF314862.gb')
        self.fasta_path = package_paths.data_path('JF314862.fasta')
        self.seq_str = 'CATCATCAACATCATCGTGCCCTGCGTGCTCATCTCCTTCGTGGCTGTGC' + \
                       'TCGTCTACTTTCTGCCTGCCAAGGGTAACGCTGGCACCAGGCGGCTGTGG' + \
                       'GACTGCCTGTGCCATAGGCGTGAAGAGGGCAGGCCATGTGGCTGGGCAGA' + \
                       'GGGAGGGAAGTGGGGGACAGCCACCGCTGGGAGACTGGCACCTGGGCCCA' + \
                       'GTGCCCGTCATTTCCCCATCACATGGGCTTGGGGACATGGAAGCCAGTCC' + \
                       'TGTGGGAGCAGACAGACACTCCCGGCTGCCGTGTCAGTCCTTAGGGCTGG' + \
                       'CTGGACTCTCTCTGCACAGCCTCCCACTGTCAGTCCCAGGACCATCCATG' + \
                       'TCCTAGGCATGTCTAGGCAGAGCCAGGCCCTTTCCAGGTGCCCTGGGACC' + \
                       'CCGTCTCACGTGTCGATCCCCTCACTCTCCACATCCTGGCAGCGGGTGGG' + \
                       'CAGAAGTGCACCGTCTCCATCAATGTCC'
        self.fasta_str = ">s1\nACGTGCTATCTATCGTATTTAG\n"
        self.small_fasta = self.getTestFile('small.fasta')
        out = open(self.small_fasta, 'w')
        out.write(self.fasta_str)
        out.close()

    def test_small_fasta(self):
        seq = read_seq(self.small_fasta, format='fasta', data_type='dna')
        self.assertIsInstance(seq, SeqRecord)
        self.assertEqual(seq.id, 's1')
        self.assertEqual(str(seq.seq), 'ACGTGCTATCTATCGTATTTAG')

    def test_fasta(self):
        seq = read_seq(self.fasta_path, format='fasta', data_type='dna')
        self.assertIsInstance(seq, SeqRecord)
        self.assertEqual(str(seq.seq), self.seq_str)

    def test_genbank(self):
        seq = read_seq(self.gb_path, format='gb', data_type='dna')
        self.assertIsInstance(seq, SeqRecord)
        self.assertEqual(str(seq.seq), self.seq_str)
        self.assertEqual(seq.name, 'JF314862')
        self.assertEqual(seq.id, 'JF314862.1')

    def test_guess_format(self):
        seq = read_seq(self.small_fasta, format=None, data_type='dna')
        self.assertIsInstance(seq, SeqRecord)
        self.assertEqual(seq.id, 's1')
        self.assertEqual(str(seq.seq), 'ACGTGCTATCTATCGTATTTAG')
        
        seq = read_seq(self.fasta_path, format=None, data_type='dna')
        self.assertIsInstance(seq, SeqRecord)
        self.assertEqual(str(seq.seq), self.seq_str)

        seq = read_seq(self.gb_path, format=None, data_type='dna')
        self.assertIsInstance(seq, SeqRecord)
        self.assertEqual(str(seq.seq), self.seq_str)
        self.assertEqual(seq.name, 'JF314862')
        self.assertEqual(seq.id, 'JF314862.1')

class GetSeqIterTestCase(SeqSiftTestCase):
    def setUp(self):
        self.gb_path = package_paths.data_path('JF314863-JF314866.gb')
        self.fasta_path = package_paths.data_path('JF314863-JF314866.fasta')
        self.names = ['JF' + str(x) for x in range(314863, 314867)]
        self.ids = [str(x) + '.1' for x in self.names]

    def test_genbank(self):
        seq_iter = get_seq_iter(self.gb_path, format='gb', data_type='dna')
        names = []
        ids = []
        for seq in seq_iter:
            self.assertIsInstance(seq, SeqRecord)
            names.append(seq.name)
            ids.append(seq.id)
        self.assertEqual(len(names), 4)
        self.assertEqual(sorted(names), sorted(self.names))
        self.assertEqual(sorted(ids),
                sorted(self.ids))

    def test_fasta(self):
        seq_iter = get_seq_iter(self.fasta_path, format='fasta',
                data_type='dna')
        names = []
        for seq in seq_iter:
            self.assertIsInstance(seq, SeqRecord)
            n = seq.name.split('|')[3]
            names.append(n)
        self.assertEqual(len(names), 4)
        self.assertEqual(sorted(names), sorted(self.ids))

    def test_genbank_fasta(self):
        si1 = get_seq_iter(self.gb_path, format='gb', data_type='dna')
        si2 = get_seq_iter(self.fasta_path, format='fasta', data_type='dna')
        self.assertSameSequenceData(si1, si2)

    def test_guess_format(self):
        si1 = get_seq_iter(self.gb_path, format=None, data_type='dna')
        si2 = get_seq_iter(self.fasta_path, format=None, data_type='dna')
        self.assertSameSequenceData(si1, si2)

class GetBufferedSeqIterTestCase(SeqSiftTestCase):
    def setUp(self):
        self.gb_path = package_paths.data_path('JF314863-JF314866.gb')
        self.seqs = get_buffered_seq_iter(self.gb_path, format='gb',
                data_type='dna')

    def test_buffer(self):
        s1 = [s for s in self.seqs]
        s2 = [s for s in self.seqs]
        self.assertSameData(s1, s2)

    def test_guess_format(self):
        sequences = get_buffered_seq_iter(self.gb_path)
        seq_list = [s for s in sequences]
        self.assertSameSequences(seq_list, self.seqs)

class GetIndexedSeqIterTestCase(SeqSiftTestCase):
    def setUp(self):
        self.gb_path = package_paths.data_path('JF314863-JF314866.gb')
        self.seqs = get_indexed_seq_iter(self.gb_path, format='gb',
                data_type='dna')

    def test_indexing(self):
        ids = []
        for id_str in self.seqs:
            self.assertIsInstance(self.seqs[id_str], SeqRecord)
            self.assertEqual(self.seqs[id_str].id, id_str)
            ids.append(id_str)
        for i in range(2):
            for id_str in ids:
                seq = self.seqs[id_str]
                self.assertIsInstance(seq, SeqRecord)
                self.assertEqual(seq.id, id_str)

class GetSeqDictTestCase(SeqSiftTestCase):
    def setUp(self):
        self.gb_path = package_paths.data_path('JF314863-JF314866.gb')
        self.seqs = get_seq_dict(self.gb_path, format='gb', data_type='dna')

    def test_dict(self):
        self.assertIsInstance(self.seqs, dict)
        for k, v in self.seqs.iteritems():
            self.assertIsInstance(v, SeqRecord)
            self.assertEqual(v.id, k)

class ConvertFormatTestCase(SeqSiftTestCase):
    def setUp(self):
        self.mkTestDir()

    def test_primates(self):
        formats = {'fasta': '.fasta', 'phylip-relaxed': '.phylip',
                'nexus': '.nexus'}
        for in_format, in_ext in formats.iteritems():
            in_file = package_paths.data_path('primates' + in_ext)
            for out_format, out_ext in formats.iteritems():
                out_file = self.getTestFile('primates' + out_ext)
                n = convert_format(in_file=in_file,
                        in_format=in_format,
                        out_file=out_file,
                        out_format=out_format,
                        data_type='dna')
                self.assertEqual(n, 12)
                in_seqs = SeqIO.parse(in_file, format=in_format,
                        alphabet=IUPAC.ambiguous_dna)
                out_seqs = SeqIO.parse(out_file, format=out_format,
                        alphabet=IUPAC.ambiguous_dna)
                self.assertSameData(in_seqs, out_seqs)

    def test_limnonectes(self):
        formats = {'fasta': '.fasta', 'phylip-relaxed': '.phylip',
                'nexus': '.nexus'}
        for in_format, in_ext in formats.iteritems():
            in_file = package_paths.data_path('limnonectes' + in_ext)
            for out_format, out_ext in formats.iteritems():
                out_file = self.getTestFile('limnonectes' + out_ext)
                n = convert_format(in_file=in_file,
                        in_format=in_format,
                        out_file=out_file,
                        out_format=out_format,
                        data_type='dna')
                self.assertEqual(n, 80)
                in_seqs = SeqIO.parse(in_file, format=in_format,
                        alphabet=IUPAC.ambiguous_dna)
                out_seqs = SeqIO.parse(out_file, format=out_format,
                        alphabet=IUPAC.ambiguous_dna)
                self.assertSameData(in_seqs, out_seqs)

    def test_caenophidia(self):
        formats = {'fasta': '.fasta', 'phylip-relaxed': '.phylip',
                'nexus': '.nexus'}
        for in_format, in_ext in formats.iteritems():
            in_file = package_paths.data_path('caenophidia' + in_ext)
            for out_format, out_ext in formats.iteritems():
                out_file = self.getTestFile('caenophidia' + out_ext)
                n = convert_format(in_file=in_file,
                        in_format=in_format,
                        out_file=out_file,
                        out_format=out_format,
                        data_type='protein')
                self.assertEqual(n, 114)
                in_seqs = SeqIO.parse(in_file, format=in_format,
                        alphabet=IUPAC.extended_protein)
                out_seqs = SeqIO.parse(out_file, format=out_format,
                        alphabet=IUPAC.extended_protein)
                self.assertSameData(in_seqs, out_seqs)

if __name__ == '__main__':
    unittest.main()
