import os
import sys
import unittest
import subprocess

from Bio import SeqIO

from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class GbFetchTestCase(SeqSiftTestCase):
    def setUp(self):
        self.gbfetch = package_paths.scripts_path("gbfetch.py")
        self.mkTestDir()
        self.out = self.getTestFile('gbfetch_test.txt')
        self.expected_fasta = package_paths.data_path('JF314863-JF314866.fasta')
        self.expected_gb = package_paths.data_path('JF314863-JF314866.gb')

    def exe_gbfetch(self, arg_list, return_code=0, stdout=None, stderr=None):
        if isinstance(arg_list, str):
            args = arg_list.split()
        else:
            args = arg_list
        cmd = ['python', self.gbfetch, '-o', self.out] + args
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
        self.exe_gbfetch(['-h'], return_code=0)

    def test_no_args(self):
        self.exe_gbfetch([], return_code=1)

    def test_fasta(self):
        self.exe_gbfetch(['-a', 'JF314863,JF314864',
                          '-g', '354698780,354698782',
                          '--format', 'fasta'])
        expected = SeqIO.parse(self.expected_fasta, format='fasta')
        results = SeqIO.parse(self.out, format='fasta')
        self.assertSameData(expected, results)

    def test_gb(self):
        self.exe_gbfetch(['-a', 'JF314863,JF314864',
                          '-g', '354698780,354698782',
                          '--format', 'gb'])
        expected = SeqIO.parse(self.expected_gb, format='gb')
        results = SeqIO.parse(self.out, format='gb')
        self.assertSameData(expected, results)
  
if __name__ == '__main__':
    unittest.main()
