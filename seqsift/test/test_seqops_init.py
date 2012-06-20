#! /usr/bin/env python

import os
import sys
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

from seqsift.seqops.__init__ import *
from seqsift.test.support import package_paths
from seqsift.test.support.extended_test_case import SeqSiftTestCase
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class CopySeqMetadataTestCase(SeqSiftTestCase):
    def setUp(self):
        self.seq = SeqIO.read(
                package_paths.data_path('JF314862.gb'),
                format='gb',
                alphabet=IUPAC.ambiguous_dna)
    
    def test_empty_seq(self):
        s = copy_seq_metadata(self.seq)
        self.assertSameMetadata(self.seq, s)
        self.assertEqual(str(s.seq), '')

    def test_basic(self):
        s = copy_seq_metadata(self.seq, 'AGCT')
        self.assertSameMetadata(self.seq, s)
        self.assertEqual(str(s.seq), 'AGCT')

if __name__ == '__main__':
    unittest.main()
