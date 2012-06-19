import os
import sys
import unittest

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

if __name__ == '__main__':
    unittest.main()
