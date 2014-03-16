#! /usr/bin/env python

import os
import sys
import argparse

from Bio import SeqIO
from seqsift.utils import FILE_FORMATS, VALID_DATA_TYPES, argparse_utils

_program_info = {
    'name': 'seq-file-slice',
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'copyright': 'Copyright (C) 2012 Jamie Oaks.',
    'license': (
        'This is free software distributed under the GNU General Public '
        'License in the hope that it will be useful, but WITHOUT ANY '
        'WARRANTY. You are free to change and redistribute it in accord with '
        'the GPL. See the GNU General Public License for more details.'),}

def main_cli():
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('input_files', metavar='INPUT-SEQ-FILE',
            nargs = '+',
            type = argparse_utils.arg_is_file,
            help = ('Input sequence file(s) to be output into files with '
                    '`-n` sequences per file.'))
    parser.add_argument('-n', '--num-seqs-per-file',
            type = int,
            required = True,
            default = 4000000,
            help = ('The maximum number of sequences to put in each output '
                    'file.'))
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
    parser.add_argument('-c', '--compress',
            action = 'store_true',
            help = 'Compress (gzip) output files.')
    parser.add_argument('-o', '--output-dir',
            type = argparse_utils.arg_is_dir,
            help = ('The directory in which all output files will be written. '
                    'The default is to use the directory of the input file.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            help = ('Prefix to use at beginning of output files. The default '
                    'is to use the first input file name.'))
    parser.add_argument('--log-frequency',
            type = argparse_utils.arg_is_nonnegative_int,
            default = 100000,
            help = ('The frequency at which to log progress. Default is to log '
                    'every 100000 sequences.'))
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

    from seqsift.utils import dataio
    from seqsift.utils.fileio import OpenFile

    ##########################################################################
    ## handle args

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

    if not args.prefix:
        args.prefix = os.path.splitext(args.input_files[0])[0]
    out_ext = FILE_FORMATS.get_ext(args.input_format,
            compressed = args.compress)

    compresslevel = None
    if args.compress:
        compresslevel = 9

    # handle sequential formats on the fly
    if FILE_FORMATS.is_sequential(args.input_format):
        seq_iter = dataio.get_seq_iter(
                file_objs = args.input_files,
                format = args.input_format,
                data_type = args.data_type)

        dataio.write_seqs_to_files(seq_iter,
                max_num_seqs_per_file = args.num_seqs_per_file,
                format = args.input_format,
                compresslevel = compresslevel,
                prefix = args.prefix,
                force = False)

    # use SeqIO for non-sequential formats
    else:
        batch_iter = dataio.get_seq_batch_iter_from_files(
                file_objs = args.input_files,
                number_per_batch = args.num_seqs_per_file,
                format = args.input_format,
                data_type = args.data_type)

        for batch_idx, seq_iter in enumerate(batch_iter):
            out_path = '{0}_{1:0>4}{2}'.format(args.prefix, batch_idx + 1,
                    out_ext)
            if os.path.exists(out_path):
                log.error('ERROR: File {0} already exists!')
                sys.exit(1)
            out = OpenFile(out_path, mode = 'w', compresslevel = compresslevel)
            SeqIO.write(seq_iter,
                    handle = out,
                    format = args.input_format)
            out.close()

if __name__ == '__main__':
    main_cli()

