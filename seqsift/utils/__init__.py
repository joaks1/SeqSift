#! /usr/bin/env python

import sys
import os
import itertools
import errno

from seqsift.utils.errors import FileExtensionError

class FileFormats(dict):
    def __init__(self):
        dict.__init__(self, {
                '.aln': 'clustal',
                '.clustal': 'clustal',
                '.fasta': 'fasta',
                '.fas': 'fasta',
                '.faa': 'fasta',
                # '.fastq': 'fastq',
                '.gb': 'genbank',
                '.gbk': 'genbank',
                '.genbank': 'genbank',
                '.nex': 'nexus',
                '.nexus': 'nexus',
                '.phy': 'phylip-relaxed',
                '.phylip': 'phylip-relaxed',
                # '.sff': 'sff',
                '.sth': 'stockholm',
                '.sto': 'stockholm',
                '.stockholm': 'stockholm',
            })

    def __str__(self):
        s = ''
        for ext in sorted(self.keys()):
            s += '\t{0:<12}: {1}\n'.format(ext, self[ext])
        return s

    def _get_supported_formats(self):
        return sorted(list(set(self.values())))

    supported_formats = property(_get_supported_formats)

    def get_format_from_extension(self, file_extension):
        try:
            return self[file_extension.lower()]
        except KeyError:
            raise FileExtensionError(
                    'Unrecognized file extension: {0!r}\n'
                    'Here are the supported options:\n{1}'.format(
                            file_extension,
                            str(self)))

    def get_format_from_path(self, file_path):
        ext = os.path.splitext(file_path)[-1]
        return self.get_format_from_extension(ext)
    
    def get_format_from_file_object(self, file_obj):
        if isinstance(file_obj, str):
            return self.get_format_from_path(file_obj)
        elif hasattr(file_obj, 'name'):
            return self.get_format_from_path(file_obj.name)
        else:
            return None

FILE_FORMATS = FileFormats()

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

VALID_DATA_TYPES = ['dna', 'rna', 'protein', 'aa']


def mkdr(path):
    """
    Creates directory `path`, but suppresses error if `path` already exists.
    """
    try:
        os.makedirs(path)
    except OSError, e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
