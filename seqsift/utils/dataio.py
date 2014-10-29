#! /usr/bin/env python

import sys
import os
import tempfile
import re
try:
    import cPickle as pickle
except ImportError:
    import pickle

from Bio.Alphabet import IUPAC
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqsift.seqops import sequtils
from seqsift.utils import (VALID_DATA_TYPES, FILE_FORMATS, functions, fileio,
        errors, alphabets)
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

class BufferedIter(object):
    def __init__(self, obj_iter = None):
        self._tmp = tempfile.TemporaryFile(mode = 'a+b')
        self._pickler = pickle.Pickler(self._tmp)
        self._unpickler = pickle.Unpickler(self._tmp)
        self.n = 0
        if obj_iter:
            self.extend(obj_iter)

    def extend(self, obj_iter):
        for obj in obj_iter:
            self.append(obj)

    def append(self, obj):
        self._pickler.dump(obj)
        self.n += 1

    def seek(self, idx = 0):
        self._tmp.seek(idx)

    def __iter__(self):
        self.seek()
        while True:
            try:
                yield self._unpickler.load()
            except EOFError:
                break

class SeqFileIter(object):
    count = 0
    def __init__(self, file_obj,
            format = None,
            data_type = 'dna',
            ambiguities = True):
        self.__class__.count += 1
        self.instance_name = '-'.join([self.__class__.__name__,
                str(self.count)])
        self.name = getattr(file_obj, 'name', self.instance_name)
        self._close = False
        self._file_obj = file_obj
        if isinstance(file_obj, str):
            self.name = file_obj
            self._file_obj = fileio.OpenFile(file_obj, 'r')
            self._close = True
        if format == None:
            format = FILE_FORMATS.get_format_from_file_object(file_obj)
        self._seqs = SeqIO.parse(self._file_obj,
                format=format,
                alphabet=get_state_alphabet(data_type, ambiguities))

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        try:
            return self._seqs.next()
        except StopIteration as e:
            if self._close:
                self._file_obj.close()
            raise e

    @classmethod
    def get_seq_iter_from_file(cls, file_obj,
            format = None,
            data_type = 'dna',
            ambiguities = True):
        """
        Returns a SeqRecord iterator from a sequence file.
        """
        return cls(file_obj,
                format = format,
                data_type = data_type,
                ambiguities = ambiguities)

    @classmethod
    def get_seq_iter(cls, file_objs,
            format = None,
            data_type = 'dna',
            ambiguities = True):
        for f in file_objs:
            seqs = cls(f,
                    format = format,
                    data_type = data_type,
                    ambiguities = ambiguities)
            for s in seqs:
                yield s


get_seq_iter_from_file = SeqFileIter.get_seq_iter_from_file
get_seq_iter = SeqFileIter.get_seq_iter

def get_buffered_seq_iter_from_file(file_obj, format=None, data_type='dna',
        ambiguities=True):
    return BufferedIter(get_seq_iter_from_file(file_obj,
            format=format,
            data_type=data_type,
            ambiguities=ambiguities))

def get_buffered_seq_iter(file_objs, format=None, data_type='dna',
        ambiguities=True):
    return BufferedIter(get_seq_iter(file_objs,
            format=format,
            data_type=data_type,
            ambiguities=ambiguities))


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
    return SeqIO.to_dict(get_seq_iter([file_obj],
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
    seqs = get_seq_iter(file_objs, format = format, data_type = data_type,
            ambiguities = ambiguities)
    return sequtils.seq_batch_iter(seqs, number_per_batch)

def write_seqs(seqs, dest = None,
        format = 'fasta',
        compresslevel = None):
    if not dest:
        dest = sys.stdout
    out, close = fileio.process_file_arg(dest, mode = 'w',
            compresslevel = compresslevel)
    for seq in seqs:
        out.write('{0}'.format(seq.format(format)))
    if close:
        out.close()

def write_seqs_to_files(seqs,
        max_num_seqs_per_file = float('inf'),
        format = 'fasta',
        compresslevel = None,
        prefix = '',
        force = False):
    compress = False
    if compresslevel:
        compress = True
    ext = FILE_FORMATS.get_ext(format, compress)
    file_idx = 0
    file_stream = None
    for seq_idx, seq in enumerate(seqs):
        if seq_idx % max_num_seqs_per_file == 0:
            if file_stream:
                file_stream.close()
            file_idx += 1
            path = '{0}_{1:0>4}{2}'.format(prefix, file_idx, ext)
            if os.path.exists(path) and (not force):
                raise errors.PathExistsError('File {0} already exists'.format(
                        path))
            file_stream = fileio.OpenFile(path, mode = 'w',
                    compresslevel = compresslevel)
        file_stream.write('{0}'.format(seq.format(format)))
    if file_stream and (not file_stream.closed):
        file_stream.close()

class LociFileIter(object):
    count = 0
    dna_alphabet = alphabets.DnaAlphabet()
    dna_symbols = ''.join(set(
            [x.upper() for x in dna_alphabet.get_valid_symbols()] +
            [x.lower() for x in dna_alphabet.get_valid_symbols()]))
    seq_pattern = re.compile(r'^>(?P<name>\S+)\s+(?P<seq>[{0}]+)$'.format(dna_symbols))
    inter_locus_pattern = re.compile(r'^//.*$')

    def __init__(self, file_obj):
        self.__class__.count += 1
        self.instance_name = '-'.join([self.__class__.__name__,
                str(self.count)])
        self.name = getattr(file_obj, 'name', self.instance_name)
        self._close = False
        self._file_obj = file_obj
        if isinstance(file_obj, str):
            self.name = file_obj
            self._file_obj = fileio.OpenFile(file_obj, 'r')
            self._close = True

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        try:
            return self._next_seq()
        except StopIteration as e:
            print "DONE"
            if self._close:
                self._file_obj.close()
            raise e

    def _next_seq(self):
        # for idx, line in enumerate(self._file_obj):
        while True:
            line = self._file_obj.next()
            m = self.seq_pattern.match(line)
            x = self.inter_locus_pattern.match(line)
            if m:
                s = Seq(m.group('seq'), alphabet = get_state_alphabet('dna',
                        ambiguities = True))
                yield SeqRecord(
                        seq = s,
                        id = m.group('name'),
                        name = m.group('name'))
            elif x:
                return
            else:
                raise Exception('unexpected format of line in loci-formatted '
                        'file {0}:\n{1}\n'.format(self.name, line.strip()))


