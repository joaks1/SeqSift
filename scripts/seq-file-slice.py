#! /usr/bin/env python

import os
import sys
import argparse

from Bio import SeqIO
from seqsift.seqops.seqfilter import column_filter, row_filter
from seqsift.utils.dataio import get_buffered_seq_iter, convert_format
from seqsift.utils import FILE_FORMATS, VALID_DATA_TYPES
from seqsift.utils.messaging import get_logger

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

_LOG = get_logger(__name__, 'INFO')

def arg_is_path(path):
    try:
        if not os.path.exists(path):
            raise
    except:
        msg = 'path {0!r} does not exist'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_file(path):
    try:
        if not is_file(path):
            raise
    except:
        msg = '{0!r} is not a file'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)

def arg_is_dir(path):
    try:
        if not is_dir(path):
            raise
    except:
        msg = '{0!r} is not a directory'.format(path)
        raise argparse.ArgumentTypeError(msg)
    return expand_path(path)


def main_cli():
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('input_files', metavar='INPUT-SEQ-FILE',
            nargs = '+',
            type = arg_is_file,
            help = ('Input sequence file(s) to be sliced up into smaller '
                    'output files.'))
    parser.add_argument('-f', '--from',
            dest = 'from_format',
            type = str,
            choices = FILE_FORMATS.supported_formats,
            help = ('The format of the input sequence file. Valid options '
                    'include: {0}. By default, the format is guessed based on '
                    'the extension of the first input file. However, if '
                    'provided, this option will always take precedence over '
                    'the file extension.'.format(
                          ', '.join(FILE_FORMATS.supported_formats))))
    parser.add_argument('-t', '--to',
            dest = 'to_format',
            type = str,
            choices = FILE_FORMATS.supported_formats,
            help = ('The desired format of the output sequence file. Valid '
                    'options include: {0}. By default, if an output file path '
                    'is provided, the format is guessed based on the extension '
                    'of this file. However, this option will always take '
                    'precedence over the file extension. Either this option or '
                    'an output file path with an extension is required; if '
                    'neither are provided the program will exit with an '
                    'error.'.format(', '.join(FILE_FORMATS.supported_formats))))
    parser.add_argument('-d', '--data-type',
            type = str,
            choices = VALID_DATA_TYPES,
            default='dna',
            help = ('The type of sequence data. The default is dna. Valid '
                    'options include: {0}.'.format(', '.join(
                            VALID_DATA_TYPES))))
    parser.add_argument('-n', '--num-seqs-per-file',
            type = int,
            required = True,
            help = ('The maximum number of sequences to put in each output '
                    'file.'))
    parser.add_argument('-o', '--output-dir',
            type = arg_is_dir,
            help = ('The directory in which all output files will be written. '
                    'The default is to use the directory of the input file.'))
    parser.add_argument('-p', '--prefix',
            action = 'store',
            type = str,
            default = '',
            help = ('Prefix to use at beginning of output files. The default '
                    'is to use the input file name.'))

    args = parser.parse_args()

    
    if not args.to_format:
        args.to_format = FILE_FORMATS.get_format_from_file_object(
                args.input_files[0])
    if not args.to_format:
        _LOG.error("Could not determine output format.\n"
                   "You must either provide the output format\n"
                   "using the '--to' option or have a recognized\n"
                   "file extension on the first input file.\n"
                   "Here are the supported file extensions:\n{0}".format(
                        str(FILE_FORMATS)))
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)

    for fp in args.input_files:
        with open(fp, 'w')
    if len(opt_dict) == 0:
        convert_format(in_file = in_file_path,
                       out_file = out_file_path,
                       in_format = in_format,
                       out_format = out_format,
                       data_type = data_type)
        sys.exit(0)

    seqs = get_buffered_seq_iter(in_file_path,
            format = in_format,
            data_type = data_type)

    if options.remove_missing_sequences:
        seqs = row_filter(seqs,
                character_list = list(options.missing_characters),
                max_frequency = options.missing_sequence_proportion)

    if options.remove_missing_columns:
        seqs = column_filter(seqs,
                character_list = list(options.missing_characters),
                max_frequency = options.missing_column_proportion)

    SeqIO.write(seqs,
                handle = out_file_path,
                format = out_format)

if __name__ == '__main__':
    main()

