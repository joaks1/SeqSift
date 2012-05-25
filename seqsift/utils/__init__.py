#! /usr/bin/env python

import itertools

DEFAULT_DNA_SIMILARITY_MATRIX = {
        ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
        ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
        ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
        ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
}

DNA_AMBIGUITY_CODES = {
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
        'N': ('A', 'C', 'G', 'T')
}
DNA_REVERSE_AMBIGUITY_CODES = {}
for k, v in DNA_AMBIGUITY_CODES.iteritems():
    for permutation in itertools.permutations(v):
        DNA_REVERSE_AMBIGUITY_CODES[permutation] = k
