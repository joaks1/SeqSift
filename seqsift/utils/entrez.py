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

LEADING_ZEROS = re.compile(r'^\s*([0]+)(\d)+\s*$')
GI = r'\s*(\d+)\s*'
GI_PATTERN = re.compile(r'^' + GI + r'$')
ACC = r'\s*([a-zA-z]{1,5})(\d{5,9})\s*'
ACC_PATTERN = re.compile(r'^' + ACC + r'$')

def parse_accession_numbers(string):
    acc_range = re.compile(r'^' + ACC + r'-\s*([a-zA-z]{0,5})(\d{5,9})\s*$')

    acc_list = string.strip().split(',')
    accs = set()
    for acc_str in acc_list:
        m1 = ACC_PATTERN    .match(acc_str.strip())
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
                accs.add(prefix1.upper() + str(i).zfill(len(num1)))
        else:
            _LOG.warning("cannot parse accession number(s) {0!r}"
                         "... skipping!".format(acc_str.strip()))
    return list(accs)

def parse_gi_numbers(string):
    gi_range = re.compile(r'^' + GI + r'-' + GI + r'$')
 
    gi_list = string.strip().split(',')
    gis = set()
    for gi_str in gi_list:
        m1 = GI_PATTERN.match(gi_str.strip())
        m2 = gi_range.match(gi_str.strip())
        if m1:
            gis.add(m1.groups()[0])
        elif m2:
            from_gi, to_gi = m2.groups()
            if int(from_gi) > int(to_gi):
                _LOG.warning("gi number range {0!r} is invalid... "
                             "skipping!".format(gi_str.strip()))
                continue
            gis.update([str(x).zfill(len(from_gi)) for x in range(int(from_gi), int(to_gi)+1)])
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

def parse_mixed_gi_list(gi_list):
    if isinstance(gi_list, str):
        ids = sorted([x.strip() for x in gi_list.split(',')], reverse=True)
    else:
        ids = sorted([x.strip() for x in gi_list], reverse=True)
    accs = [x.strip() for x in ids if ACC_PATTERN.match(x)]
    gis = [x.strip() for x in ids if GI_PATTERN.match(x)]
    return list(set(gis)), list(set(accs))

def get_pure_gi_numbers(gi_list):
    gis, accs = parse_mixed_gi_list(gi_list)
    new_gis = []
    if accs:
         _LOG.info('searching for gi numbers for accessions {0}'.format(accs))
    for a in accs:
        s = Entrez.read(Entrez.esearch(db='nuccore', term=a+"[accession]"))
        if ('IdList' in s) and (len(s['IdList']) == 1):
            _LOG.info('found gi number ({0}) for {1}'.format(
                    s['IdList'][0], a))
            new_gis.append(s['IdList'][0])
        elif ('IdList' in s) and (len(s['IdList']) > 1):
            _LOG.warn('found multiple gi numbers ({0}) for {1}'.format(
                    s['IdList'], a))
            new_gis.extend(s['IdList'])
        elif ('IdList' in s) and (len(s['IdList']) < 1):
            _LOG.warn('could not find gi number for {0}'.format(a))
        else:
            _LOG.warn('could not find gi number for {0}'.format(a))
    return sorted(list(set(gis + new_gis)))

def entrez_post(gi_list, db):
    ids = get_pure_gi_numbers(gi_list)
    post = Entrez.read(Entrez.epost(db=db, id=','.join(ids)))
    return post["WebEnv"], post["QueryKey"]
        
def get_gb_handle(gi_list, db, rettype, retmode='text', tmp_file=False):
    _LOG.info("fetching genbank ids {0}".format(gi_list))
    if isinstance(gi_list, str):
        ids = [x.strip() for x in gi_list.split(',')]
    else:
        ids = gi_list
    if len(ids) > 10:
        webenv, query_key = entrez_post(ids, db=db)
        h = Entrez.efetch(db=db, webenv=webenv, query_key=query_key,
                rettype=rettype, retmode=retmode)
    else:
        h = Entrez.efetch(db=db, id=','.join(ids), rettype=rettype,
                retmode=retmode)
    if tmp_file:
        return dataio.get_tmp_handle(h, rewind=True)
    else:
        return h

def fetch_gb_seqs(gi_list, data_type,
        parse_function=dataio.get_seq_iter_from_file):
    db = get_entrez_database(data_type)
    file_obj = get_gb_handle(gi_list, db=db, rettype='gb', retmode='text')
    return parse_function(file_obj, format='gb', data_type=data_type)
