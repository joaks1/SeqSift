#! /usr/bin/env python

import sys
import os
import random

from seqsift.utils.errors import FileExtensionError
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class FileFormats(dict):
    def __init__(self):
        dict.__init__(self, {
                '.aln': 'clustal',
                '.clustal': 'clustal',
                '.fasta': 'fasta',
                '.fas': 'fasta',
                '.faa': 'fasta',
                '.fastq': 'fastq',
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
        _extensions = {
                'clustal': '.aln',
                'fasta': '.fasta',
                'fastq': '.fastq',
                'genbank': '.gb',
                'nexus': '.nex',
                'phylip-relaxed': '.phy',
                'phylip': '.phy',
                'stockholm': '.sto'}

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
        if self.has_gzip_ext(file_path):
            file_path = os.path.splitext(file_path)[0]
        ext = os.path.splitext(file_path)[-1]
        return self.get_format_from_extension(ext)
    
    def get_format_from_file_object(self, file_obj):
        if isinstance(file_obj, str):
            return self.get_format_from_path(file_obj)
        elif hasattr(file_obj, 'name'):
            return self.get_format_from_path(file_obj.name)
        else:
            return None

    def has_gzip_ext(self, file_path):
        return (os.path.splitext(file_path)[-1].lower() == '.gz')

    def get_ext(self, file_format, compressed = False):
        ext = ''
        try:
            ext = self._extensions[file_format]
        except KeyError as e:
            _LOG.error('unsupported file format {0}'.format(file_format))
            raise e
        if compressed:
            return ext + '.gz'
        return ext

FILE_FORMATS = FileFormats()

DEFAULT_DNA_SIMILARITY_MATRIX = {
        ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
        ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
        ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
        ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
}

VALID_DATA_TYPES = ['dna', 'rna', 'protein', 'aa']

GLOBAL_RNG = random.Random()

