#! /usr/bin/env python

import os
import unittest
import tempfile
import cPickle
import re

from seqsift.dataio import BufferedIter
from seqsift.test.support import package_paths
from seqsift.utils import mkdr

class SeqSiftTestCase(unittest.TestCase):
    test_paths = []
    test_dir = None

    def mkTestDir(self):
        self.test_dir = package_paths.output_path()
        mkdr(self.test_dir)

    def getTestFile(self, filename):
        tp = package_paths.output_path(filename)
        self.test_paths.append(tp)
        return tp

    def getTestStream(self, filename):
        tp = self.getTestFile(filename)
        return open(tp, 'w')

    def getBufferedIter(self, sequences):
        tmp = tempfile.TemporaryFile()
        pickler = cPickle.Pickler(tmp)
        for seq in sequences:
            pickler.dump(seq)

        def seq_iter(sk=0):
            tmp.seek(sk)
            unpickler = cPickle.Unpickler(tmp)
            while True:
                try:
                    yield unpickler.load()
                except EOFError:
                    break

        yield seq_iter

    def remove_gaps(self, seq_str):
        return re.sub(r'[-]', '', seq_str)

    def assertSameIDs(self, sequences1, sequences2):
        self.assertEqual(
                sorted([s.id for s in sequences1]),
                sorted([s.id for s in sequences2]))

    def assertSameNames(self, sequences1, sequences2):
        self.assertEqual(
                sorted([s.name for s in sequences1]),
                sorted([s.name for s in sequences2]))

    def assertSameDescriptions(self, sequences1, sequences2):
        self.assertEqual(
                sorted([s.description for s in sequences1]),
                sorted([s.description for s in sequences2]))

    def assertSameSequenceData(self, sequences1, sequences2):
        seqs1 = [self.remove_gaps(str(s.seq)) for s in sequences1]
        seqs2 = [self.remove_gaps(str(s.seq)) for s in sequences2]
        self.assertEqual(sorted(seqs1), sorted(seqs2))

    def assertSameSequences(self, sequences1, sequences2, aligned=False):
        if isinstance(sequences1, BufferedIter):
            s1 = sequences1
        else:
            s1 = BufferedIter(sequences1)
        if isinstance(sequences2, BufferedIter):
            s2 = sequences2
        else:
            s2 = BufferedIter(sequences2)
        for seq1 in s1.iter():
            try:
                seq2 = next(s for s in s2.iter() if s.id == seq1.id)
            except StopIterations:
                fail("failed assertSameSequences:\n"
                     "Sequence id {0!r} of {1!r} not found in {2!r}\n".format(
                        seq1.id, sequences1, sequences2))
            if aligned:
                self.assertEqual(str(seq1.seq), str(seq2.seq))
            else:
                self.assertEqual(
                        self.remove_gaps(str(seq1.seq)),
                        self.remove_gaps(str(seq2.seq)))

    def assertSameData(self, sequences1, sequences2, aligned=False):
        self.assertSameIDs(sequences1, sequences2)
        self.assertSameNames(sequences1, sequences2)
        self.assertSameDescriptions(sequences1, sequences2)
        self.assertSameSequences(sequences1, sequences2, aligned=aligned)

    def tearDown(self):
        for p in self.test_paths:
            if os.path.exists(p):
                os.remove(p)
        if self.test_dir and os.path.exists(self.test_dir):
            os.removedirs(self.test_dir)
