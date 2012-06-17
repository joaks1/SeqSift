#! /usr/bin/env python

import sys
import os
import tempfile
import cPickle
import contextlib
import warnings

from Bio.Alphabet import IUPAC
from Bio import SeqIO, AlignIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqsift.utils import VALID_DATA_TYPES
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)
warnings.filterwarnings(action="ignore", category=UserWarning,
        module=r'.*Entrez.*')

class BufferedIter(object):
    def __init__(self, obj_iter):
        self._tmp = tempfile.TemporaryFile()
        self._pickler = cPickle.Pickler(self._tmp)
        self._unpickler = cPickle.Unpickler(self._tmp)
        for obj in obj_iter:
            self._pickler.dump(obj)

    def iter(self):
        self._tmp.seek(0)
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

def read_seq(file_obj, format, data_type, ambiguities=True):
    """
    Returns a single SeqRecord from a file containing exactly one sequence
    record.
    """
    _LOG.info("reading sequence from {0!r}.".format(file_obj))
    return SeqIO.read(file_obj,
            format=format,
            alphabet=get_state_alphabet(data_type, ambiguities))

def get_seq_iter(file_obj, format, data_type, ambiguities=True):
    """
    Returns a SeqRecord iterator from a sequence file.
    """
    _LOG.info("parsing SeqRecord iterator from {0!r}".format(file_obj))
    return SeqIO.parse(file_obj,
            format=format,
            alphabet=get_state_alphabet(data_type, ambiguities))

def get_buffered_seq_iter(file_obj, format, data_type, ambiguities=True):
    return BufferedIter(get_seq_iter(file_obj,
            format=format,
            data_type=data_type,
            ambiguities=ambiguities))

def get_indexed_seq_iter(file_path, format, data_type, key_function=None,
        ambiguities=True):
    """
    Returns an indexed SeqRecord iterator from a sequence file.
    Only supports sequential file formats (e.g., fasta and genbank).
    The iterator acts like a dict, but only parses a sequence from
    a file when needed. For large sequence files, this is memory-efficient
    alternative to reading all the sequences into a dict.
    """
    _LOG.info("parsing indexed SeqRecord iterator from {0!r}.".format(
            file_path))
    return SeqIO.index(file_path,
            format=format,
            alphabet=get_state_alphabet(data_type, ambiguities),
            key_function=key_function)

def get_seq_dict(file_obj, format, data_type, ambiguities=True):
    """
    Returns a dict of SeqRecords from a sequenc file.
    This loads all the sequences in the file into memory. This is
    efficient for small sequence files, but may cause memory issues
    with large files.
    """
    return SeqIO.to_dict(get_seq_iter(file_obj,
            format=format,
            data_type=data_type,
            ambiguities=ambiguities))

def convert_format(in_file, in_format, out_file, out_format, data_type,
        ambiguities=True):
    _LOG.info("converting {in_format}-formatted file {in_file!r} to "
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

def get_entrez_database(data_type):
    dt = data_type.lower()
    if dt == 'dna' or dt == 'rna':
        return 'nuccore'
    elif dt == 'protein' or dt == 'aa':
        return 'protein'
    else:
        raise ValueError(
                "'{0!r}' is not a valid data type. Options:\n\t{1}".format(
                        data_type, ", ".join(VALID_DATA_TYPES)))

def get_gb_handle(gi_list, db, rettype, retmode='text', tmp_file=False):
    _LOG.info("fetching genbank ids {0}".format(gi_list))
    if isinstance(gi_list, str):
        ids = gi_list
    else:
        ids = ",".join(gi_list)
    if len(gi_list) > 10:
        post = Entrez.read(Entrez.epost(db=db, id=ids))
        webenv = post["WebEnv"]
        query_key = post["QueryKey"]
        h = Entrez.efetch(db=db, webenv=webenv, query_key=query_key,
                rettype=rettype, retmode=retmode)
    else:
        h = Entrez.efetch(db=db, id=ids, rettype=rettype, retmode=retmode)
    if tmp_file:
        return get_tmp_handle(h, rewind=True)
    else:
        return h

def fetch_gb_seqs(gi_list, data_type, parse_function=get_seq_iter):
    db = get_entrez_database(data_type)
    file_obj = get_gb_handle(gi_list, db=db, rettype='gb', retmode='text')
    return parse_function(file_obj, format='gb', data_type=data_type)
