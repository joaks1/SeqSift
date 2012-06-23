#! /usr/bin/env python

import os
import sys
from optparse import OptionParser

from seqsift.utils.entrez import (get_gb_handle, parse_gi_numbers, 
        parse_accession_numbers)
from seqsift.utils.messaging import get_logger

_program_info = {
    'name': 'gbfetch',
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'copyright': 'Copyright (C) 2012 Jamie Oaks.',
    'license': (
        'This is free software distributed under the GNU General Public '
        'License in the hope that it will be useful, but WITHOUT ANY '
        'WARRANTY. You are free to change and redistribute it in accord with '
        'the GPL. See the GNU General Public License for more details.'),}

_LOG = get_logger(__name__, 'INFO')

_FORMATS = {'fas': 'fasta',
            'fasta': 'fasta',
            'gb': 'gb',
            'genbank': 'gb',}

def main():
    description = '{name} {version}'.format(**_program_info)
    usage = ("\n  %prog [options]")
    parser = OptionParser(usage=usage, description=description,
                          version=_program_info['version'],
                          add_help_option=True)
    parser.add_option("-a", "--accessions", dest="accessions", 
        type="string",
        help=("GenBank accession numbers. "
              "E.g., -g JF314862,JF314864-314866 -OR-"
              " -g 'JF314862, JF314864-314866'"))
    parser.add_option("-g", "--gi_numbers", dest="gi_numbers", 
        type="string",
        help=("GenBank GI numbers. E.g., -g 354698774,354698776-354698779 -OR-"
              " -g '354698774, 354698776-354698779'"))
    parser.add_option("--format", dest="format", type="string",
        default='gb',
        help=("Desired file format. Valid options are 'fasta' and 'gb'. "
              "Default is 'gb'."))
    parser.add_option("-d", "--db", dest="entrez_database", type="string",
        default='nuccore',
        help=("The Entrez database to search. Default is 'nuccore'. "
              "Must be a valid Entrez E-utility Database name: "
              "http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2."
              "chapter2_table1/?report=objectonly"))
    parser.add_option("-o", "--output", dest="output_filepath", type="string",
        help="Path to output file. Default is standard output")
    (options, args) = parser.parse_args()

    if not options.accessions and not options.gi_numbers:
        _LOG.error("You must provide GI or accession numbers")
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)
    if not options.format.lower() in _FORMATS.keys():
        _LOG.error("Format {0} is not valid".format(options.format))
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)
    format = _FORMATS[options.format.lower()]
    if options.output_filepath:
        out = open(options.output_filepath, 'w')
    else:
        out = sys.stdout
    
    gi_list = []
    if options.gi_numbers:
        gi_list = parse_gi_numbers(options.gi_numbers)
    if options.accessions:
        gi_list += parse_accession_numbers(options.accessions)

    h = get_gb_handle(gi_list=gi_list,
                      db=options.entrez_database,
                      rettype=format,
                      retmode='text',
                      tmp_file=False)
    for line in h:
        out.write(line)
    h.close()
    if options.output_filepath:
        out.close()

if __name__ == '__main__':
    main()
