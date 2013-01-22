#! /usr/bin/env python

import os
import unittest
import tempfile
import cPickle
import re

from seqsift.utils.dataio import BufferedIter
from seqsift.test.support import package_paths
from seqsift.utils import mkdr
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class SeqSiftTestCase(unittest.TestCase):
    test_paths = []
    test_dir = None

    def mkTestDir(self):
        self.test_dir = package_paths.output_path()
        mkdr(self.test_dir)

    def getTestFile(self, filename):
        tp = package_paths.output_path(filename)
        self.appendTestFile(tp)
        return tp
    
    def appendTestFile(self, filepath):
        _LOG.debug("appending {0!r} to test paths".format(filepath))
        self.test_paths.append(filepath)

    def getTestStream(self, filename):
        tp = self.getTestFile(filename)
        return open(tp, 'w'), tp

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

    def assertSameSequenceData(self, sequences1, sequences2, aligned=False):
        if aligned:
            seqs1 = [str(s.seq) for s in sequences1]
            seqs2 = [str(s.seq) for s in sequences2]
        else:
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
        for seq1 in s1:
            try:
                seq2 = next(s for s in s2 if s.id == seq1.id)
            except StopIteration:
                self.fail("failed assertSameSequences:\n"
                     "Sequence id {0!r} of {1!r} not found in {2!r}\n".format(
                        seq1.id, sequences1, sequences2))
            if aligned:
                self.assertEqual(str(seq1.seq), str(seq2.seq))
            else:
                self.assertEqual(
                        self.remove_gaps(str(seq1.seq)),
                        self.remove_gaps(str(seq2.seq)))

    def assertSameData(self, sequences1, sequences2, aligned=False,
            include_descriptions=False):
        if isinstance(sequences1, BufferedIter):
            s1 = sequences1
        else:
            s1 = BufferedIter(sequences1)
        if isinstance(sequences2, BufferedIter):
            s2 = sequences2
        else:
            s2 = BufferedIter(sequences2)
        self.assertSameIDs(s1, s2)
        self.assertSameNames(s1, s2)
        if include_descriptions:
            self.assertSameDescriptions(s1, s2)
        self.assertSameSequences(s1, s2, aligned=aligned)

    def assertSameMetadata(self, seq1, seq2):
        self.assertEqual(seq1.id, seq2.id)
        self.assertEqual(seq1.name, seq2.name)
        self.assertEqual(seq1.description, seq2.description)
        self.assertEqual(seq1.letter_annotations, seq2.letter_annotations)
        self.assertEqual(seq1.annotations, seq2.annotations)
        self.assertEqual(seq1.features, seq2.features)
        self.assertEqual(seq1.dbxrefs, seq2.dbxrefs)

    def tearDown(self):
        for p in self.test_paths:
            if os.path.exists(p):
                _LOG.debug("removing test path {0!r}".format(p))
                os.remove(p)
        if self.test_dir and os.path.exists(self.test_dir):
            os.removedirs(self.test_dir)
