import os
import sys
import unittest
import subprocess
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)
sys.path = [package_paths.scripts_path()] + sys.path
seqdigest = __import__('seqdigest')

class testParseGiNumbersTestCase(unittest.TestCase):
    def test_empty_string(self):
        gis = seqdigest.parse_gi_numbers('')
        self.assertEqual(gis, [])
    
    def test_invalid_gis(self):
        gis = seqdigest.parse_gi_numbers('AA123, 221!#3, 123A34-123C344')
        self.assertEqual(gis, [])

    def test_invalid_range(self):
        gis = seqdigest.parse_gi_numbers('1234-1223')
        self.assertEqual(gis, [])
    
    def test_basic(self):
        gis = seqdigest.parse_gi_numbers('1212, 1234, 2300-2312')
        self.assertEqual(sorted(gis), [1212, 1234] + range(2300,2313))


class SeqDigestTestCase(SeqSiftTestCase):
    DATA_LINE = re.compile(r'^(\d+)\t(\d+)\n$')
    def setUp(self):
        self.seqdigest = package_paths.scripts_path("seqdigest.py")
        self.mkTestDir()        

    def parse_result_file(self, file_path):
        len_dist = {}
        for line in open(file_path, 'rU'):
            m = self.DATA_LINE.match(line)
            if m:
                l, f = m.groups()
                len_dist[int(l)] = int(f)
        return len_dist

    def exe_seqdigest(self, arg_list, return_code=0, stdout=None, stderr=None):
        if isinstance(arg_list, str):
            args = arg_list.split()
        else:
            args = arg_list
        cmd = ['python', self.seqdigest, '-o', self.test_dir] + args
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
        self.exe_seqdigest(['-h'], return_code=0)

    def test_no_args(self):
        self.exe_seqdigest([], return_code=1)

    def test_basic(self):
        expected = {
                'JF314863': {
                        13: 2,
                        65: 1,
                        105: 1,
                        117: 1,
                        175: 1,},
                'JF314864': {
                        13: 1,
                        65: 1,
                        115: 1,
                        117: 1,
                        180: 1,},
                'JF314865': {
                        13: 1,
                        75: 1,
                        96: 1,
                        103: 1,
                        178: 1,},
                'JF314866': {
                        13: 2,
                        65: 1,
                        105: 1,
                        117: 1,
                        175: 1,},
                'combined' : {
                        13: 6,
                        65: 3,
                        75: 1,
                        96: 1,
                        103: 1,
                        105: 2,
                        115: 1,
                        117: 3,
                        175: 2,
                        178: 1,
                        180: 1,}}
        rs = 'TAG'
        cs = '3'
        self.exe_seqdigest(['-s', rs,
                            '-c', cs,
                            '-g', '354698776,354698778',
                            package_paths.data_path('JF314865-JF314866.gb')])
        results = {}
        for k in expected.iterkeys():
            result_file_path = os.path.join(self.test_dir,
                    ".".join([k, 'txt']))
            self.appendTestFile(result_file_path)
            results[k] = self.parse_result_file(result_file_path)
        self.assertEqual(expected, results)

    def test_extra_length(self):
        expected = {
                'JF314863': {
                        23: 2,
                        75: 1,
                        115: 1,
                        127: 1,
                        185: 1,},
                'JF314864': {
                        23: 1,
                        75: 1,
                        125: 1,
                        127: 1,
                        190: 1,},
                'JF314865': {
                        23: 1,
                        85: 1,
                        106: 1,
                        113: 1,
                        188: 1,},
                'JF314866': {
                        23: 2,
                        75: 1,
                        115: 1,
                        127: 1,
                        185: 1,},
                'combined' : {
                        23: 6,
                        75: 3,
                        85: 1,
                        106: 1,
                        113: 1,
                        115: 2,
                        125: 1,
                        127: 3,
                        185: 2,
                        188: 1,
                        190: 1,}}
        rs = 'TAG'
        cs = '3'
        self.exe_seqdigest(['-s', rs,
                            '-c', cs,
                            '-g', '354698776,354698778',
                            '-x', '10',
                            package_paths.data_path('JF314865-JF314866.gb')])
        results = {}
        for k in expected.iterkeys():
            result_file_path = os.path.join(self.test_dir,
                    ".".join([k, 'txt']))
            self.appendTestFile(result_file_path)
            results[k] = self.parse_result_file(result_file_path)
        self.assertEqual(expected, results)
        
if __name__ == '__main__':
    unittest.main()
