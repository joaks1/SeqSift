#! /usr/bin/env python

"""
A CLI tool for summarizing sequence length information.

The program will generate a brief summary (min, max, mean, variance) of the
length of sequences found in one or more files.
"""

import os
import sys
import argparse

from Bio import SeqIO
from seqsift.utils import FILE_FORMATS, VALID_DATA_TYPES, argparse_utils

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2014 Jamie Oaks.',
    'license': (
        'This is free software distributed under the GNU General Public '
        'License in the hope that it will be useful, but WITHOUT ANY '
        'WARRANTY. You are free to change and redistribute it in accord with '
        'the GPL. See the GNU General Public License for more details.'),}

def write_summary(name, summary_dict, stream = sys.stdout):
    stream.write('[{0}]\n'.format(name))
    stream.write('\tn = {0}\n'.format(summary_dict.n))
    stream.write('\tmin length = {0}\n'.format(summary_dict.minimum))
    stream.write('\tmax length = {0}\n'.format(summary_dict.maximum))
    stream.write('\tmean length = {0}\n'.format(summary_dict.mean))
    stream.write('\tlength variance = {0}\n'.format(summary_dict.variance))
    stream.write('\tlength standard deviation = {0}\n'.format(summary_dict.std_deviation))

def main_cli():
    description = '{name} {version}\n\n{description}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description,
            formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_files', metavar='INPUT-SEQ-FILE',
            nargs = '+',
            type = argparse_utils.arg_is_file,
            help = ('Input sequence file(s) to be output into files with '
                    '`-n` sequences per file.'))
    parser.add_argument('--format',
            dest = 'input_format',
            type = str,
            choices = FILE_FORMATS.supported_formats,
            help = ('The format of the input sequence file(s). Valid options '
                    'include: {0}. By default, the format is guessed based on '
                    'the extensions of input file(s). However, if '
                    'provided, this option will always take precedence over '
                    'the file extension.'.format(
                          ', '.join(FILE_FORMATS.supported_formats))))
    parser.add_argument('-d', '--data-type',
            type = str,
            choices = VALID_DATA_TYPES,
            default='dna',
            help = ('The type of sequence data. The default is dna. Valid '
                    'options include: {0}.'.format(', '.join(
                            VALID_DATA_TYPES))))
    parser.add_argument('--quiet',
            action = 'store_true',
            help = 'Run without verbose messaging.')
    parser.add_argument('--debug',
            action = 'store_true',
            help = 'Run in debugging mode.')

    args = parser.parse_args()

    ##########################################################################
    ## set up logging

    from seqsift.utils.messaging import get_logger, LOGGING_LEVEL_ENV_VAR

    os.environ[LOGGING_LEVEL_ENV_VAR] = "INFO"
    if args.quiet:
        os.environ[LOGGING_LEVEL_ENV_VAR] = "WARNING"
    if args.debug:
        os.environ[LOGGING_LEVEL_ENV_VAR] = "DEBUG"
    log = get_logger(name = __name__)

    ##########################################################################
    ## package imports

    from seqsift.seqops import seqstats

    ##########################################################################
    ## handle args

    if not args.input_format:
        args.input_format = None

    summaries = seqstats.get_seq_summaries_from_files(args.input_files,
            format = args.input_format,
            data_type =  args.data_type)
    global_summary = summaries.pop('global')

    keys = sorted(summaries.keys())
    for k in keys:
        write_summary(k, summaries[k])
    write_summary('overall', global_summary)


if __name__ == '__main__':
    main_cli()

