#! /usr/bin/env python

import sys
import os
import tempfile
import cPickle
from itertools import islice, chain

from Bio.Alphabet import IUPAC
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqsift.seqops import sequtils
from seqsift.utils import VALID_DATA_TYPES, FILE_FORMATS, functions
from seqsift.utils.messaging import get_logger
from seqsift.utils.errors import FileExtensionError

_LOG = get_logger(__name__)

class BufferedIter(object):
    def __init__(self, obj_iter = None):
        self._tmp = tempfile.TemporaryFile(mode = 'a+b')
        self._pickler = cPickle.Pickler(self._tmp)
        self._unpickler = cPickle.Unpickler(self._tmp)
        if obj_iter:
            self.extend(obj_iter)

    def extend(self, obj_iter):
        for obj in obj_iter:
            self.append(obj)

    def append(self, obj):
        self._pickler.dump(obj)

    def seek(self, idx = 0):
        self._tmp.seek(idx)

    def __iter__(self):
        self.seek()
        while True:
            try:
                yield self._unpickler.load()
            except EOFError:
                break

def get_tmp_handle(file_obj=None, rewind=True):
    tmp = tempfile.TemporaryFile()
    if not file_obj:
        return tmp
    elif isinstance(file_obj, str):
        f = open(file_obj, 'rU')
    else:
        f = file_obj
    for line in f:
        tmp.write(line)
    if rewind:
        tmp.seek(0)
    return tmp

def get_state_alphabet(data_type, ambiguities=True):
    if data_type.lower() == 'dna':
        if ambiguities:
            return IUPAC.ambiguous_dna
        else:
            return IUPAC.unambiguous_dna
    elif data_type.lower() == 'rna':
        if ambiguities:
            return IUPAC.ambiguous_rna
        else:
            return IUPAC.unambiguous_rna
    elif data_type.lower() == 'protein' or data_type.lower() == 'aa':
        if ambiguities:
            return IUPAC.extended_protein
        else:
            return IUPAC.protein
    else:
        raise ValueError(
                "'{0!r}' is not a valid data type. Options:\n\t{1}".format(
                        data_type, ", ".join(VALID_DATA_TYPES)))

def read_seq(file_obj, format=None, data_type='dna', ambiguities=True):
    """
    Returns a single SeqRecord from a file containing exactly one sequence
    record.
    """
    if format == None:
        format = FILE_FORMATS.get_format_from_file_object(file_obj)
    _LOG.debug("reading sequence from {0!r}.".format(file_obj))
    return SeqIO.read(file_obj,
            format=format,
            alphabet=get_state_alphabet(data_type, ambiguities))

def get_seq_iter(file_obj, format=None, data_type='dna', ambiguities=True):
    """
    Returns a SeqRecord iterator from a sequence file.
    """
    if format == None:
        format = FILE_FORMATS.get_format_from_file_object(file_obj)
    _LOG.debug("parsing SeqRecord iterator from {0!r}".format(file_obj))
    return SeqIO.parse(file_obj,
            format=format,
            alphabet=get_state_alphabet(data_type, ambiguities))

def seq_iter(file_objs, format = None, data_type = 'dna', ambiguities = True):
    for f in file_objs:
        seqs = get_seq_iter(f, format=format, data_type=data_type,
                ambiguities=ambiguities)
        for s in seqs:
            yield s

def get_buffered_seq_iter(file_obj, format=None, data_type='dna',
        ambiguities=True):
    if format == None:
        format = FILE_FORMATS.get_format_from_file_object(file_obj)
    return BufferedIter(get_seq_iter(file_obj,
            format=format,
            data_type=data_type,
            ambiguities=ambiguities))

def buffered_seq_iter(file_objs, format=None, data_type='dna',
        ambiguities=True):
    return BufferedIter(seq_iter(file_objs,
            format=format,
            data_type=data_type,
            ambiguities=ambiguities))

def get_indexed_seq_iter(file_path, format=None, data_type='dna',
        key_function=None,
        ambiguities=True):
    """
    Returns an indexed SeqRecord iterator from a sequence file.  Only supports
    sequential file formats (e.g., fasta and genbank).  The iterator acts like
    a dict, but only parses a sequence from a file when needed. For large
    sequence files, this is memory-efficient alternative to reading all the
    sequences into a dict.
    """
    if format == None:
        format = FILE_FORMATS.get_format_from_file_object(file_path)
    _LOG.debug("parsing indexed SeqRecord iterator from {0!r}.".format(
            file_path))
    return SeqIO.index(file_path,
            format=format,
            alphabet=get_state_alphabet(data_type, ambiguities),
            key_function=key_function)

def get_seq_dict(file_obj, format=None, data_type='dna', ambiguities=True):
    """
    Returns a dict of SeqRecords from a sequence file.  This loads all the
    sequences in the file into memory. This is efficient for small sequence
    files, but may cause memory issues with large files.
    """
    if format == None:
        format = FILE_FORMATS.get_format_from_file_object(file_obj)
    return SeqIO.to_dict(get_seq_iter(file_obj,
            format=format,
            data_type=data_type,
            ambiguities=ambiguities))

def convert_format(in_file, out_file,
        in_format=None,
        out_format=None,
        data_type='dna',
        ambiguities=True):
    if in_format == None:
        in_format = FILE_FORMATS.get_format_from_file_object(in_file)
    if out_format == None:
        out_format = FILE_FORMATS.get_format_from_file_object(out_file)
    _LOG.debug("converting {in_format}-formatted file {in_file!r} to "
              "{out_format}-formatted file {out_file!r}.".format(
                    in_file=in_file,
                    in_format=in_format,
                    out_file=out_file,
                    out_format=out_format))
    nseqs = SeqIO.convert(
            in_file=in_file,
            in_format=in_format,
            out_file=out_file,
            out_format=out_format,
            alphabet=get_state_alphabet(data_type, ambiguities))
    return nseqs

def get_seq_batch_iter_from_files(file_objs,
        number_per_batch,
        format = None,
        data_type = 'dna',
        ambiguities = True):
    seqs = seq_iter(file_objs, format = format, data_type = data_type,
            ambiguities = ambiguities)
    return sequtils.seq_batch_iter(seqs, number_per_batch)

