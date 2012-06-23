#! /usr/bin/env python

import sys
import os
import re
import warnings

from Bio import Entrez

from seqsift.utils import dataio
from seqsift.utils import VALID_DATA_TYPES
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)
warnings.filterwarnings(action="ignore", category=UserWarning,
        module=r'.*Entrez.*')

def parse_accession_numbers(string):
    acc = r'\s*([a-zA-z]{1,5})(\d{5,9})\s*'
    acc_pattern = re.compile(r'^' + acc + r'$')
    acc_range = re.compile(r'^' + acc + r'-\s*([a-zA-z]{0,5})(\d{5,9})\s*$')

    acc_list = string.strip().split(',')
    accs = set()
    for acc_str in acc_list:
        m1 = acc_pattern.match(acc_str.strip())
        m2 = acc_range.match(acc_str.strip())
        if m1:
            prefix, num = m1.groups()
            accs.add(prefix.upper() + num)
        elif m2:
            prefix1, num1, prefix2, num2 = m2.groups()
            if prefix2 and not (prefix1.upper() == prefix2.upper()):
                _LOG.warning("invalid accession range {0!r}; "
                             "cannot range across prefixes"
                             "...skipping!".format(acc_str.strip()))
                continue
            if int(num1) > int(num2):
                _LOG.warning("accession range {0!r} is invalid... "
                             "skipping!".format(acc_str.strip()))
                continue
            for i in range(int(num1), int(num2)+1):
                accs.add(prefix1.upper() + str(i))
        else:
            _LOG.warning("cannot parse accession number(s) {0!r}"
                         "... skipping!".format(acc_str.strip()))
    return list(accs)

def parse_gi_numbers(string):
    gi = r'\s*(\d+)\s*'
    gi_pattern = re.compile(r'^' + gi + r'$')
    gi_range = re.compile(r'^' + gi + r'-' + gi + r'$')
 
    gi_list = string.strip().split(',')
    gis = set()
    for gi_str in gi_list:
        m1 = gi_pattern.match(gi_str.strip())
        m2 = gi_range.match(gi_str.strip())
        if m1:
            gis.add(m1.groups()[0])
        elif m2:
            from_gi, to_gi = m2.groups()
            if int(from_gi) > int(to_gi):
                _LOG.warning("gi number range {0!r} is invalid... "
                             "skipping!".format(gi_str.strip()))
                continue
            gis.update([str(x) for x in range(int(from_gi), int(to_gi)+1)])
        else:
            _LOG.warning("cannot parse gi number(s) {0!r}... skipping!".format(
                    gi_str.strip()))
    return list(gis)

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
        return dataio.get_tmp_handle(h, rewind=True)
    else:
        return h

def fetch_gb_seqs(gi_list, data_type, parse_function=dataio.get_seq_iter):
    db = get_entrez_database(data_type)
    file_obj = get_gb_handle(gi_list, db=db, rettype='gb', retmode='text')
    return parse_function(file_obj, format='gb', data_type=data_type)
