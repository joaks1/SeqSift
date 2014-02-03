#! /usr/bin/env python

import os
import sys
from optparse import OptionParser, OptionGroup

from Bio import SeqIO
from seqsift.seqops import seqmod
from seqsift.seqops.seqfilter import column_filter, row_filter
from seqsift.utils.dataio import get_buffered_seq_iter, convert_format
from seqsift.utils import FILE_FORMATS, VALID_DATA_TYPES
from seqsift.utils.messaging import get_logger

_program_info = {
    'name': 'seqaid',
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'copyright': 'Copyright (C) 2012 Jamie Oaks.',
    'license': (
        'This is free software distributed under the GNU General Public '
        'License in the hope that it will be useful, but WITHOUT ANY '
        'WARRANTY. You are free to change and redistribute it in accord with '
        'the GPL. See the GNU General Public License for more details.'),}

_LOG = get_logger(__name__, 'INFO')

def main():
    description = '{name} {version}'.format(**_program_info)
    usage = ("\n  %prog [options] <SEQ_INPUT_FILE> [<SEQ_OUTPUT_FILE>]")
    parser = OptionParser(usage=usage, description=description,
                          version=_program_info['version'],
                          add_help_option=True)
    format_opts = OptionGroup(parser, 'Format Options',
            'These options designate file formats and data type.')
    format_opts.add_option('-f', '--from', dest='from_format', type='string',
            help=('The format of the input sequence file. Valid options '
                  'include: {0}. By default, the format is guessed based on '
                  'the extension of the input file. However, if provided, '
                  'this option will always take precedence over the file '
                  'extension.'.format(
                        ', '.join(FILE_FORMATS.supported_formats))))
    format_opts.add_option('-t', '--to', dest='to_format', type='string',
            help=('The desired format of the output sequence file. Valid '
                  'options include: {0}. By default, if an output file path '
                  'is provided, the format is guessed based on the extension '
                  'of this file. However, this option will always take '
                  'precedence over the file extension. Either this option or '
                  'an output file path with an extension is required; if '
                  'neither are provided the program will exit with an '
                  'error.'.format(', '.join(FILE_FORMATS.supported_formats))))
    format_opts.add_option('-d', '--data-type', dest='data_type', type='string',
            default='dna',
            help=('The type of sequence data. The default is dna. Valid '
                  'options include: {0}.'.format(', '.join(VALID_DATA_TYPES))))
    parser.add_option_group(format_opts)

    filter_opts = OptionGroup(parser, 'Filter Options',
            'These options allow filtering of data by columns or sequences.')
    filter_opts.add_option('--remove-missing-columns',
            dest='remove_missing_columns',
            default=False,
            action='store_true',
            help=("Remove aligned columns with missing data. Characters to be "
                  "considered missing can be specified with the "
                  "--missing-characters option; the default is '?-'. "
                  "The proportion of rows that must contain these characters "
                  "for a row to be removed can be specified with the "
                  "--missing-column-proportion option; the default is 1.0. "
                  "Note, this option is only relevant to aligned sequences, "
                  "and will result in an error if the input sequences are not "
                  "aligned."))
    filter_opts.add_option('--missing-column-proportion',
            dest='missing_column_proportion',
            type='float',
            default=1.0,
            help=('The proportion of rows that must contain '
                  '--missing-characters for a column to be removed. '
                  'This option is only relevant in combination with the '
                  '--remove-missing-columns option.'))
    filter_opts.add_option('--remove-missing-sequences',
            dest='remove_missing_sequences',
            default=False,
            action = 'store_true',
            help=("Remove sequences with missing data. Characters to be "
                  "considered missing can be specified with the "
                  "--missing-characters option; the default is '?-'. "
                  "The proportion of the sites that must contain these "
                  "characters for a sequence to be removed can be specified "
                  "with the --missing-sequence-proportion option; the default "
                  "is 1.0."))
    filter_opts.add_option('--missing-sequence-proportion',
            dest='missing_sequence_proportion',
            type='float',
            default=1.0,
            help=('The proportion of sites that must contain '
                  '--missing-characters for a sequence to be removed. '
                  'This option is only relevant in combination with the '
                  '--remove-missing-sequences option.'))
    filter_opts.add_option('--missing-characters', dest='missing_characters',
            type='str',
            default='?-',
            help=("Characters to be considered missing and be used in "
                  "evaluating columns/sequences to remove with the "
                  "--remove-missing-columns and --remove-missing-sequences "
                  "options. The default is '?-'."))
    parser.add_option_group(filter_opts)

    rev_comp_opts = OptionGroup(parser, 'Reverse Complement Options',
            'These options are for reverse complementing sequences.')
    rev_comp_opts.add_option('--rev-comp',
            dest='rev_comp',
            default = False,
            action = 'store_true',
            help=("Reverse complement all sequences. This option overrides "
                  "all other reverse-complement options."))
    rev_comp_opts.add_option('--fix-rev-comp-by',
            dest='fix_rev_comp_by',
            type = 'choice',
            choices = ['first', 'read'],
            help=("Try to correct reverse complement errors. "
                  "Options include 'first' and 'read'. If 'first' is "
                  "specified, sequences are returned in their orientation "
                  "that minimizes distance from the first sequence. "
                  "If 'read' is used, sequences are returned in their "
                  "orientation that has the longest read frame "
                  "(see 'Translation Options' for controlling translation "
                  "of reading frames)."))
    parser.add_option_group(rev_comp_opts)

    translation_opts = OptionGroup(parser, 'Translation Options',
            ('These options control translation from nucleotide to amino acid '
             'sequences.'))
    translation_opts.add_option('--table',
            type = 'choice',
            choices = list(range(1, 7)) + list(range(9, 17)) + list(range(21, 26)),
            default = 1,
            help = ('The translation table to use for any options associated '
                    'with translating nucleotide sequences to amino acids. '
                    'Option should be the integer that corresponds to the '
                    'desired translation table according to NCBI '
                    '(http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). '
                    'The default is 1 (the "standard" code).'))
    translation_opts.add_option('--allow-partial',
        default = False,
        action = 'store_true',
        help = ('Allow partial reading frames at the beginning (no start '
                'codon) and end (no stop codon) of sequences.'))
    translation_opts.add_option('--read-after-stop',
        default = False,
        action = 'store_true',
        help = ('A new reading frame begins immediately after a stop codon. '
                'The default is to start reading frame at next start codon '
                'after a stop codon. This option might be useful for exons.'))
    parser.add_option_group(translation_opts)
    (options, args) = parser.parse_args()
    
    if len(args) == 1:
        in_file_path = args[0]
        out_file_path = sys.stdout
    elif len(args) == 2:
        in_file_path = args[0]
        out_file_path = args[1]
    elif len(args) > 2:
        _LOG.error("Too many arguments. Expecting at most 2 arguments:\n"
                   "The path to the input file (required), and the path to\n"
                   "output file (optional; defaults to standard output).")
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)
    elif len(args) < 1:
        _LOG.error("Too few arguments. Expecting at least 1 argument:\n"
                   "the path to the input file.")
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)

    opt_dict = options.__dict__

    if options.from_format:
        in_format = opt_dict.pop('from_format')
    else:
        in_format = FILE_FORMATS.get_format_from_file_object(in_file_path)
    if not in_format:
        _LOG.error("Could not determine format of input file.\n"
                   "You must either provide the format of the input file\n"
                   "using the '--from-format' option or have a recognized\n"
                   "file extension on the input file. Here are the supported\n"
                   "file extensions:\n{0}".format(str(FILE_FORMATS)))
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)

    if options.to_format:
        out_format = opt_dict.pop('to_format')
    else:
        out_format = FILE_FORMATS.get_format_from_file_object(out_file_path)
    if not out_format:
        _LOG.error("Could not determine format of output file.\n"
                   "You must either provide the format of the output file\n"
                   "using the '--to-format' option or have a recognized\n"
                   "file extension on the output file. Here are the supported\n"
                   "file extensions:\n{0}".format(str(FILE_FORMATS)))
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)

    data_type = opt_dict.pop('data_type')

    if len(opt_dict) == 0:
        convert_format(in_file = in_file_path,
                       out_file = out_file_path,
                       in_format = in_format,
                       out_format = out_format,
                       data_type = data_type)
        sys.exit(0)

    if ((options.rev_comp or options.fix_rev_comp_by) and
            (data_type.lower() not in ['dna', 'rna'])):
        _LOG.error("You have selected an option for reverse complementing\n"
                   "sequences but the data type is not DNA or RNA.")
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)

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

    if options.rev_comp:
        _LOG.info('Reverse complementing all sequences...')
        seqs = seqmod.reverse_complement(seqs)
    elif options.fix_rev_comp_by == 'first':
        _LOG.info('Reverse complementing to match first sequence...')
        seqs = seqmod.reverse_complement_to_first_seq(seqs,
                per_site = False,
                aligned = False,
                ignore_gaps = True,
                alphabet = None,
                aligner_tools = ['muscle', 'mafft'],
                log_frequency = 100)
    elif options.fix_rev_comp_by == 'read':
        _LOG.info('Reverse complementing to longest reading frame...')
        seqs = seqmod.reverse_complement_to_longest_reading_frame(seqs,
                gap_characters=['-'],
                table = options.table,
                allow_partial = options.allow_partial,
                require_start_after_stop = (not options.read_after_stop),
                log_frequency = 100)

    SeqIO.write(seqs,
                handle = out_file_path,
                format = out_format)

if __name__ == '__main__':
    main()

