#! /usr/bin/env python

import os
import sys
import unittest

from seqsift.utils import alphabets, iteritems
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class TestResidueAlphabet(unittest.TestCase):
    def test_gap_init_error(self):
        self.assertRaises(ValueError, alphabets.ResidueAlphabet,
                'ACGT',
                '-')

    def test_ambiguity_init_error(self):
        self.assertRaises(ValueError, alphabets.ResidueAlphabet,
                'AB',
                {'N': ('A', 'C')})

    def test_simple(self):
        a = alphabets.ResidueAlphabet(states = 'ACGT')
        self.assertEqual(a.states, ('A', 'C', 'G', 'T'))
        self.assertEqual(a.states, a.standardize_states(a.states))
        self.assertEqual(a.gap, None)
        self.assertEqual(a.missing, None)
        self.assertEqual(a.ambiguity_codes, {})
        self.assertEqual(a.reverse_ambiguity_codes, {})
        self.assertEqual(a.residue_ambiguity_codes, {})
        self.assertEqual(a.all_residue_codes, {'A': ('A',), 'C': ('C',),
                'G': ('G',), 'T': ('T',)})

    def test_gap_error(self):
        a = alphabets.ResidueAlphabet(states = 'ACGT')
        self.assertRaises(ValueError, a._set_gap_state, '-')
        a.states = 'ACGT-'
        a.gap = '-'
        self.assertEqual(a.gap, '-')

    def test_ambiguity_error(self):
        a = alphabets.ResidueAlphabet(states = 'ACGT-')
        self.assertRaises(ValueError, a._set_ambiguity_codes,
                {'?': ('A', 'C', 'G', 'U', '-')})
        a.ambiguity_codes = {'?': 'ACGT-'}
        self.assertEqual(a.ambiguity_codes, {'?': tuple(sorted('ACGT-'))})

class TestDnaAlphabet(unittest.TestCase):
    def setUp(self):
        self.states = tuple(sorted('ACGT-'))
        self.ambiguity_codes = {
                'R': ('A', 'G'),
                'Y': ('C', 'T'),
                'K': ('G', 'T'),
                'M': ('A', 'C'),
                'S': ('C', 'G'),
                'W': ('A', 'T'),
                'V': ('A', 'C', 'G'),
                'H': ('A', 'C', 'T'),
                'D': ('A', 'G', 'T'),
                'B': ('C', 'G', 'T'),
                'N': ('A', 'C', 'G', 'T'),
                '?': self.states,
                }
        self.reverse_ambiguity_codes = {}
        for k, v in iteritems(self.ambiguity_codes):
            self.reverse_ambiguity_codes[v] = k
        self.residue_ambiguity_codes = {
                'R': ('A', 'G'),
                'Y': ('C', 'T'),
                'K': ('G', 'T'),
                'M': ('A', 'C'),
                'S': ('C', 'G'),
                'W': ('A', 'T'),
                'V': ('A', 'C', 'G'),
                'H': ('A', 'C', 'T'),
                'D': ('A', 'G', 'T'),
                'B': ('C', 'G', 'T'),
                'N': ('A', 'C', 'G', 'T'),
                }
        self.all_residue_codes = {
                'A': ('A',),
                'C': ('C',),
                'G': ('G',),
                'T': ('T',),
                }
        for k, v in iteritems(self.residue_ambiguity_codes):
            self.all_residue_codes[k] = v

    def test_dna_alphabet(self):
        a = alphabets.DnaAlphabet()
        self.assertEqual(self.states, a.states)
        self.assertEqual(self.ambiguity_codes, a.ambiguity_codes)
        self.assertEqual(self.residue_ambiguity_codes, a.residue_ambiguity_codes)
        self.assertEqual(self.all_residue_codes, a.all_residue_codes)
        self.assertEqual(a.gap, '-')
        self.assertEqual(a.missing, '?')
        self.assertEqual(a.get_symbol('AG'), 'R')
        self.assertEqual(a.get_symbol('AGCT-'), '?')
        self.assertEqual(a.get_symbol('AGCT'), 'N')




if __name__ == '__main__':
    unittest.main()

