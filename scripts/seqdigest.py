#! /usr/bin/env python

import os
import sys
import re
import itertools
import logging
from optparse import OptionParser

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from seqsift.digest import Fragment, RecognitionSeq, DigestSummary
from seqsift.utils.dataio import fetch_gb_seqs
from seqsift.utils.messaging import get_logger

_program_info = {
    'name': 'SeqDigest',
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'copyright': 'Copyright (C) 2012 Jamie Oaks.',
    'license': (
        'This is free software distributed under the GNU General Public '
        'License in the hope that it will be useful, but WITHOUT ANY '
        'WARRANTY. You are free to change and redistribute it in accord with '
        'the GPL. See the GNU General Public License for more details.'),}

_LOG = get_logger(__name__)

GI_PATTERN = re.compile(r'^\s*(\d+)\s*$')
GI_RANGE = re.compile(r'^\s*(\d+)\s*-\s*(\d+)\s*$')
EXTENSIONS = {'fas': 'fasta',
              'fasta': 'fasta',
              'gb': 'gb',
              'genbank': 'gb',}

def parse_gi_numbers(string):
    gi_list = string.strip().split(',')
    gis = set()
    for gi_str in gi_list:
        m1 = GI_PATTERN.match(gi_str.strip())
        m2 = GI_RANGE.match(gi_str.strip())
        if m1:
            gis.add(int(m1.groups()[0]))
        elif m2:
            from_gi, to_gi = m2.groups()
            if int(from_gi) > int(to_gi):
                _LOG.warning("gi number range {0!r} is invalid... "
                             "skipping!".format(gi_str.strip()))
            gis.update(list(range(int(from_gi), int(to_gi)+1)))
        else:
            _LOG.warning("cannot parse gi number(s) {0!r}... skipping!".format(
                    gi_str.strip()))
    return list(gis)

def digest_seq(recognition_seq, seq_record):
    _LOG.info("Digesting seq {0} with recognition seq {1}...".format(
            seq_record.id, str(recognition_seq.seq)))
    return DigestSummary(recognition_seq, seq_record)
        
         
def main():
    description = '{name} {version}'.format(**_program_info)
    usage = "%prog [options] [<GENBANK_FILE1> <GENBANK_FILE2> ...]"
    parser = OptionParser(usage=usage, description=description,
                          version=_program_info['version'],
                          add_help_option=True)
    parser.add_option("-v", "--verbose", dest="verbose", default=False, 
        action="store_true",
        help="Verbose output.")
    parser.add_option("-d", "--debugging", dest="debugging", default=False, 
        action="store_true",
        help="Run in debugging mode.")
    parser.add_option("-s", "--recognition_seq", dest="recognition_seq",
        type="string",
        help="Recognition sequence of restriction enzyme.")
    parser.add_option("-g", "--gi_numbers", dest="gi_numbers", 
        type="string",
        help=("GenBank GI numbers. E.g., -g 354698774,354698776-354698779 -OR-"
              " -g '354698774, 354698776-354698779'"))
    parser.add_option("--format", dest="format", type="string",
        help=("Format of sequence files. Valid options are 'fasta' and 'gb', "
              "for fasta- and genbank-formatted files, respectively. The "
              "default is to guess the format based on your file extensions; "
              ".fas or .fasta translate to 'fasta', whereas .gb or .genbank "
              "translate to 'gb'. This option overrides the default behavior, "
              "in which case, all files must be of the format provided with "
              "this flag."))
    parser.add_option("--min_length", dest="min_length", type="int", default=0,
        help="Minimum fragment length to include in count.")
    parser.add_option("--max_length", dest="max_length", type="int",
        help="Maximum fragment length to include in count.")
    (options, args) = parser.parse_args()

    if options.debugging:
        _LOG.setLevel(logging.DEBUG)
    elif options.verbose and not options.debugging:
        _LOG.setLevel(logging.INFO)
    else:
        _LOG.setLevel(logging.WARNING)
    
    if not options.recognition_seq:
        _LOG.error("You must provide a recognition sequence")
        sys.stderr.write(parser.print_help())
        sys.exit(1)
    if not options.gi_numbers and len(args) < 1:
        _LOG.error("You must provide a sequence to digest, either via the "
                   "gi number option or sequence file arguments")
        sys.stderr.write(parser.print_help())
        sys.exit(1)
    format=None
    if options.format:
        format = EXTENSIONS[options.format.lower()]

    rseq = RecognitionSeq(str(options.recognition_seq))
    gi_list = []
    if options.gi_numbers:
        gi_list = parse_gi_numbers(options.gi_numbers)
    
    digests = {}
    for gi in gi_list:
        _LOG.info("Downloading gi {0}...".format(gi))
        seq_iter = fetch_gb_seqs(str(gi), data_type='dna')
        for seq in seq_iter:
            digests[seq.id] = digest_seq(rseq, seq)

    for file_path in args:
        try:
            f = open(file_path, 'rU')
        except IOError, e:
            _LOG.error("Could not open file {0}... skipping!".format(
                    file_path))
            continue
        if not options.format:
            format = EXTENSIONS[file_path.split('.')[-1].strip().lower()]
        seq_iter = get_seq_iter(file_obj=f,
                                format=format,
                                data_type='dna',
                                ambiguities=True)
        for seq in seq_iter:
            if seq.id in digests.keys():
                _LOG.warning(
                        "Sequence {0} already digested... skipping!".format(
                        seq.id))
                continue
            digests[seq.id] = digest_seq(rseq, seq)
        
if __name__ == '__main__':
    main()
