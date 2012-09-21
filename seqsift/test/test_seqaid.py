#! /usr/bin/env python

import os
import sys
import unittest
import subprocess
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils import FILE_FORMATS
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class SeqAidTestCase(SeqSiftTestCase):
    def setUp(self):
        self.seqaid = package_paths.scripts_path("seqaid.py")
        self.mkTestDir()        
        self.simple_alignment = [
                SeqRecord(Seq('ACGT?', alphabet=IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ACGT-', alphabet=IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ACGT?', alphabet=IUPAC.ambiguous_dna), id='3'),
                SeqRecord(Seq('ACGT-', alphabet=IUPAC.ambiguous_dna), id='4'),
                SeqRecord(Seq('ACGT?', alphabet=IUPAC.ambiguous_dna), id='5')]
        stream, self.simple_alignment_path = self.getTestStream(
                'simple_alignment.fasta')
        SeqIO.write(self.simple_alignment, stream, format='fasta')
        stream.close()

    def exe_seqaid(self, arg_list, return_code=0, stdout=None, stderr=None):
        if isinstance(arg_list, str):
            args = arg_list.split()
        else:
            args = arg_list
        cmd = ['python', self.seqaid] + args
        p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        o, e = p.communicate()
        r = p.wait()
        if r != return_code:
            _LOG.error("return code {0} did not match {1}".format(r,
                    return_code))
            _LOG.error("here is the stdout:\n{0}".format(o))
            _LOG.error("here is the stderr:\n{0}".format(e))
        self.assertEquals(r, return_code)
        if stdout != None:
            self.assertEqual(o, stdout)
        if stderr != None:
            self.asserEqual(e, stderr)
    
    def test_help(self):
        self.exe_seqaid(['-h'], return_code=0)

    def test_no_args(self):
        self.exe_seqaid([], return_code=1)

    def test_three_args(self):
        self.exe_seqaid(['one', 'two', 'three'], return_code=1)

    def test_format_conversion(self):
        for in_ext, in_format in FILE_FORMATS.iteritems():
            test_stream, test_path = self.getTestStream('simple' + in_ext)
            SeqIO.write(self.simple_alignment, test_stream, format=in_format)
            test_stream.close()
            for out_ext, out_format in FILE_FORMATS.iteritems():
                if out_ext == in_ext:
                    continue
                out_path = self.getTestFile(test_path.replace(in_ext, out_ext))
                self.exe_seqaid([test_path, out_path])
                seqs_in = SeqIO.parse(test_path, format=in_format,
                        alphabet=IUPAC.ambiguous_dna)
                seqs_out = SeqIO.parse(out_path, format=out_format,
                        alphabet=IUPAC.ambiguous_dna)
                self.assertSameSequences(seqs_in, seqs_out, aligned=True)

    def test_format_conversion_protein(self):
        for filename in ['caenophidia.fasta', 'caenophidia.phylip',
                'caenophidia.nexus']:
            in_ext = os.path.splitext(filename)[-1]
            if in_ext == '.phylip':
                in_format = 'phylip-relaxed'
            else:
                in_format = in_ext.replace('.', '')
            in_path = package_paths.data_path(filename)
            for out_ext, out_format in FILE_FORMATS.iteritems():
                if out_ext == in_ext:
                    continue
                if out_format == 'genbank':
                    continue
                out_path = self.getTestFile(filename.replace(in_ext, out_ext))
                _LOG.info('converting {0} to {1}'.format(filename, out_ext))
                self.exe_seqaid([in_path, out_path, '-d', 'aa'])
                seqs_in = SeqIO.parse(in_path, format=in_format,
                        alphabet=IUPAC.extended_protein)
                seqs_out = SeqIO.parse(out_path, format=out_format,
                        alphabet=IUPAC.extended_protein)
                self.assertSameSequenceData(seqs_in, seqs_out, aligned=True)

    def test_remove_missing_column(self):
        out_path = self.getTestFile('test.fasta')
        self.exe_seqaid([
                self.simple_alignment_path,
                out_path,
                '--remove-missing-columns',
                "--missing-characters='?-'",
                '--missing-column-proportion=1.0'])
        results = SeqIO.parse(out_path, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        expected = [
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='3'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='4'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='5')]
        self.assertSameSequences(results, expected, aligned=True)

        self.exe_seqaid([
                self.simple_alignment_path,
                out_path,
                '--remove-missing-columns',
                "--missing-characters='a'",
                '--missing-column-proportion=1.0'])
        results = SeqIO.parse(out_path, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        expected = [
                SeqRecord(Seq('CGT?', alphabet=IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('CGT-', alphabet=IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('CGT?', alphabet=IUPAC.ambiguous_dna), id='3'),
                SeqRecord(Seq('CGT-', alphabet=IUPAC.ambiguous_dna), id='4'),
                SeqRecord(Seq('CGT?', alphabet=IUPAC.ambiguous_dna), id='5')]
        self.assertSameSequences(results, expected, aligned=True)

    def test_remove_partial_missing_column(self):
        out_path = self.getTestFile('test.fasta')
        self.exe_seqaid([
                self.simple_alignment_path,
                out_path,
                '--remove-missing-columns',
                "--missing-characters='?'",
                '--missing-column-proportion=0.6'])
        results = SeqIO.parse(out_path, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        expected = [
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='3'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='4'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='5')]
        self.assertSameSequences(results, expected, aligned=True)

        self.exe_seqaid([
                self.simple_alignment_path,
                out_path,
                '--remove-missing-columns',
                "--missing-characters='-'",
                '--missing-column-proportion=0.4'])
        results = SeqIO.parse(out_path, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        expected = [
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='1'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='2'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='3'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='4'),
                SeqRecord(Seq('ACGT', alphabet=IUPAC.ambiguous_dna), id='5')]
        self.assertSameSequences(results, expected, aligned=True)

        self.exe_seqaid([
                self.simple_alignment_path,
                out_path,
                '--remove-missing-columns',
                "--missing-characters='?'",
                '--missing-column-proportion=0.61'])
        results = SeqIO.parse(out_path, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameSequences(results, self.simple_alignment, aligned=True)

        self.exe_seqaid([
                self.simple_alignment_path,
                out_path,
                '--remove-missing-columns',
                "--missing-characters='-'",
                '--missing-column-proportion=0.41'])
        results = SeqIO.parse(out_path, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameSequences(results, self.simple_alignment, aligned=True)

    def test_remove_row(self):
        out_path = self.getTestFile('test.fasta')
        self.exe_seqaid([
                self.simple_alignment_path,
                out_path,
                '--remove-missing-sequences',
                "--missing-characters='-'",
                '--missing-sequence-proportion=0.2'])
        results = SeqIO.parse(out_path, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameSequences(results, self.simple_alignment[::2],
                aligned=True)

        self.exe_seqaid([
                self.simple_alignment_path,
                out_path,
                '--remove-missing-sequences',
                "--missing-characters='?'",
                '--missing-sequence-proportion=0.2'])
        results = SeqIO.parse(out_path, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameSequences(results, self.simple_alignment[1::2],
                aligned=True)

        self.exe_seqaid([
                self.simple_alignment_path,
                out_path,
                '--remove-missing-sequences',
                "--missing-characters='?'",
                '--missing-sequence-proportion=0.201'])
        results = SeqIO.parse(out_path, format='fasta',
                alphabet=IUPAC.ambiguous_dna)
        self.assertSameSequences(results, self.simple_alignment,
                aligned=True)

if __name__ == '__main__':
    unittest.main()
