#! /usr/bin/env python

"""
A CLI tool for randomly sub-sampling sequences.

The program will randomly sample sequences, without replacement, from all of
the sequences provided in one or more sequence files. The random sample of
sequences is output to standard output (or redirected file) in the same format
as the first input file.
"""

import os
import sys
import random
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

def main_cli():
    description = '{name} {version}\n\n{description}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description,
            formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_files', metavar='INPUT-SEQ-FILE',
            nargs = '+',
            type = argparse_utils.arg_is_file,
            help = ('Input sequence file(s) from which to randomly sub-sample '
                    'sequences (without replacement).'))
    parser.add_argument('-n', '--num-samples',
            type = int,
            required = True,
            help = ('The number of sequences to randomly sample.'))
    parser.add_argument('--format',
            dest = 'input_format',
            type = str,
            choices = FILE_FORMATS.supported_formats,
            help = ('The format of the input sequence file(s). Valid options '
                    'include: {0}. By default, the format is guessed based on '
                    'the extension of the first input file. However, if '
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
    parser.add_argument('--seed',
            action = 'store',
            type = int,
            help = ('Random number seed.'))
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

    from seqsift.utils import dataio, GLOBAL_RNG, functions

    ##########################################################################
    ## handle args

    ## set seed if randomly sampling sequences
    if not args.seed:
        args.seed = random.randint(1, 999999999)
    GLOBAL_RNG.seed(args.seed)
    log.warning('Seed: {0}'.format(args.seed))

    if not args.input_format:
        args.input_format = FILE_FORMATS.get_format_from_file_object(
                args.input_files[0])
    if not args.input_format:
        log.error("Could not determine input format.\n"
                   "You must either provide the input format\n"
                   "using the '--from' option or have a recognizable\n"
                   "file extension on the first input file.\n"
                   "Here are the supported file extensions:\n{0}".format(
                        str(FILE_FORMATS)))
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)

    seqs = dataio.get_seq_iter(args.input_files,
            format = args.input_format,
            data_type = args.data_type)
    samples = functions.sample_iter(iterable = seqs,
            sample_size = args.num_samples)

    SeqIO.write(samples,
            handle = sys.stdout,
            format = args.input_format)

if __name__ == '__main__':
    main_cli()

