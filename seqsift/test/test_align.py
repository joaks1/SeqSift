#! /usr/bin/env python

import unittest

from seqsift.align import *

class MatchScoreTestCase(unittest.TestCase):
    def setUp(self):
        self.similarity_matrix = {
            ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
        }
        self.asymmetric_sim_matrix = {
            ('A', 'A'): 10, ('A', 'G'):  -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -2, ('G', 'G'):   7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -6, ('C', 'G'): -10, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -8, ('T', 'G'):  -6, ('T', 'C'):  0, ('T', 'T'):  8,
        }

    def test_base_error(self):
        self.assertRaises(ValueError, match_score,
                'A', 'Z', self.similarity_matrix)
        self.assertRaises(ValueError, match_score,
                'Z', 'T', self.similarity_matrix)
        self.assertRaises(ValueError, match_score,
                'X', 'Z', self.similarity_matrix)

    def test_unambiguous_base_score(self):
        for i in ['A', 'C', 'G', 'T']:
            for j in ['A', 'C', 'G', 'T']:
                score = match_score(i, j, self.similarity_matrix)
                self.assertEqual(score, self.similarity_matrix[i, j])

    def test_ambiguous_base_score(self):
        score = match_score('N', 'N', self.similarity_matrix)
        self.assertAlmostEqual(score, 0.125, 6)
        score = match_score('N', 'A', self.similarity_matrix)
        self.assertAlmostEqual(score, 0.5, 6)

    def test_asymmetric_base_score(self):
        for i in ['A', 'C', 'G', 'T']:
            for j in ['A', 'C', 'G', 'T']:
                score = match_score(i, j, self.asymmetric_sim_matrix)
                self.assertEqual(score, self.asymmetric_sim_matrix[i, j])
        score = match_score('N', 'A', self.asymmetric_sim_matrix)
        self.assertAlmostEqual(score, -1.5, 6)
        score = match_score('A', 'N', self.asymmetric_sim_matrix)
        self.assertAlmostEqual(score, 0.5, 6)

class CalculateFMatrixTestCase(unittest.TestCase):
    def setUp(self):
        self.similarity_matrix = {
            ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
        }
        self.gap_cost = -5

    def test_small_no_gaps(self):
        f = calculate_F_matrix('ACGT', 'ACGT',
                similarity_matrix = self.similarity_matrix,
                gap_cost = self.gap_cost)
        correct = [[  0, -5, -10, -15, -20],
                   [ -5, 10,   5,   0,  -5],
                   [-10,  5,  19,  14,   9],
                   [-15,  0,  14,  26,  21],
                   [-20, -5,   9,  21,  34]]
        self.assertEqual(f, correct)

    def test_small_gap(self):
        f = calculate_F_matrix('ACCGT', 'ACGT',
                similarity_matrix = self.similarity_matrix,
                gap_cost = self.gap_cost)
                    #      A    C    G    T
        correct = [[  0,  -5, -10, -15, -20],
                   [ -5,  10,   5,   0,  -5],# A
                   [-10,   5,  19,  14,   9],# C
                   [-15,   0,  14,  14,  14],# C
                   [-20,  -5,   9,  21,  16],# G
                   [-25, -10,   4,  16,  29]]# T
        self.assertEqual(f, correct)

    def test_small_ambiguous_base(self):
        f = calculate_F_matrix('ACGT', 'ANGT',
                similarity_matrix = self.similarity_matrix,
                gap_cost = self.gap_cost)
                    #     A       C      G       T
        correct = [[  0, -5,    -10,   -15,    -20],
                   [ -5, 10,      5,     0,     -5],# A
                   [-10,  5,  10.25,  5.25,   0.25],# N
                   [-15,  0,   5.25, 17.25,  12.25],# G
                   [-20, -5,   0.25, 12.25,  25.25]]# T
        for i in range(len(f)):
            for j in range(len(f[i])):
                self.assertAlmostEqual(f[i][j], correct[i][j], 6)

class TraceMaxScoreTestCase(unittest.TestCase):
    def setUp(self):
        self.similarity_matrix = {
            ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
        }
        self.gap_cost = -5

    def test_seq_matrix_confict(self):
        fmatrix = [[  0, -5, -10, -15, -20],
                   [ -5, 10,   5,   0,  -5],
                   [-10,  5,  19,  14,   9],
                   [-15,  0,  14,  26,  21],
                   [-20, -5,   9,  21,  34]]
        self.assertRaises(Exception, trace_max_score,
                fmatrix=fmatrix, 
                seq1='AAAA',
                seq2='AAAA',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
    def test_small_no_gaps(self):
        fmatrix = [[  0, -5, -10, -15, -20],
                   [ -5, 10,   5,   0,  -5],
                   [-10,  5,  19,  14,   9],
                   [-15,  0,  14,  26,  21],
                   [-20, -5,   9,  21,  34]]
        a = trace_max_score(fmatrix=fmatrix, 
                seq1='ACGT',
                seq2='ACGT',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
        self.assertEqual(('ACGT', 'ACGT'), a)

    def test_small_gap(self):
                    #      A    C    G    T
        fmatrix = [[  0,  -5, -10, -15, -20],
                   [ -5,  10,   5,   0,  -5],# A
                   [-10,   5,  19,  14,   9],# C
                   [-15,   0,  14,  14,  14],# C
                   [-20,  -5,   9,  21,  16],# G
                   [-25, -10,   4,  16,  29]]# T
        a = trace_max_score(fmatrix=fmatrix, 
                seq1='ACCGT',
                seq2='ACGT',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
        self.assertEqual(('ACCGT', 'A-CGT'), a)
                    #      A    C    C    G    T
        fmatrix = [[  0,  -5, -10, -15, -20, -25],
                   [ -5,  10,   5,   0,  -5, -10],# A
                   [-10,   5,  19,  14,   9,   4],# C
                   [-15,   0,  14,  14,  21,  16],# G
                   [-20,  -5,   9,  14,  16,  29]]# T
        a = trace_max_score(fmatrix=fmatrix, 
                seq1='ACGT',
                seq2='ACCGT',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
        self.assertEqual(('A-CGT', 'ACCGT'), a)

    def test_small_ambiguous_base(self):
                    #     A       C      G       T
        fmatrix = [[  0, -5,    -10,   -15,    -20],
                   [ -5, 10,      5,     0,     -5],# A
                   [-10,  5,  10.25,  5.25,   0.25],# N
                   [-15,  0,   5.25, 17.25,  12.25],# G
                   [-20, -5,   0.25, 12.25,  25.25]]# T
        a = trace_max_score(fmatrix=fmatrix, 
                seq1='ACGT',
                seq2='ANGT',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
        self.assertEqual(('ACGT', 'ANGT'), a)

class GlobalAlignTestCase(unittest.TestCase):
    def setUp(self):
        self.similarity_matrix = {
            ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
        }
        self.gap_cost = -5

    def test_small_no_gap(self):
        a = global_align(seq1='ACGT', seq2='ACGT',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
        self.assertEqual(('ACGT', 'ACGT'), a)

    def test_small_gap(self):
        a = global_align(seq1='ACGT', seq2='ACCGT',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
        self.assertEqual(('A-CGT', 'ACCGT'), a)
        a = global_align(seq1='ACCGT', seq2='ACGT',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
        self.assertEqual(('ACCGT', 'A-CGT'), a)

    def test_small_ambiguous_base(self):
        a = global_align(seq1='ACGT', seq2='ANGT',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
        self.assertEqual(('ACGT', 'ANGT'), a)
        a = global_align(seq1='ANGT', seq2='ACGT',
                similarity_matrix=self.similarity_matrix,
                gap_cost=self.gap_cost)
        self.assertEqual(('ANGT', 'ACGT'), a)
