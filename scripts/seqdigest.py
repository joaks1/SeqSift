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

from seqsift.digest import Fragment, RecognitionSeq
from seqsift.utils.dataio import *
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
        help="Recognition sequence of restriction enzyme")
    parser.add_option("-g", "--gi_numbers", dest="gi_numbers", 
        type="string",
        help=("GenBank GI numbers. E.g., -g 354698774,354698776-354698779 -OR-"
              " -g '354698774, 354698776-354698779'"))
    (options, args) = parser.parse_args()

    if options.debugging:
        _LOG.setLevel(logging.DEBUG)
    elif options.verbose and not options.debugging:
        _LOG.setLevel(logging.INFO)
    else:
        _LOG.setLevel(logging.WARNING)

    gi_list = None
    if options.gi_numbers:
        gi_list = parse_gi_numbers(options.gi_numbers)

if __name__ == '__main__':
    main()
