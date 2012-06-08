import sys
import os

from seqsift.dataio import get_gb_handle, get_seq_iter
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

def digest_from_gb(gi_list, recognition_seq, tmp_files=False):
    if isinstance(gi_list, str):
        ids = (i.strip() for i in gi_list.split(','))
    else:
        ids = gi_list
    for gi in ids:
        file_obj = get_gb_handle(str(gi), db='nuccore', rettype='gb',
                retmode='text', tmp_file=tmp_files)

def digest_from_file(file_list, recognition_seq):
    for f in file_list:
        if isinstance(f, str):
            file_obj = open(f, 'rU')
        else:
            file_obj = f
        seqs = get_seq_iter(f, format='gb', data_type='dna', ambiguities=True)
        digest_seq(seq.seq, recognition_seq)

def digest_seq(seq, recognition_seq):
    last_cut_index = -len(recognition_seq)
    for i in range(0, len(seq) - len(recognition_seq) + 1):
        if str(seq[i: i + len(recognition_seq)]) == str(recognition_seq):
            print i
            if (i - last_cut_index) < len(recognition_seq):
                print "overlap!"
            last_cut_index = i

class Digest(object):
    def __init__(self):
        pass
