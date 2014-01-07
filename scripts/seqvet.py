#! /usr/bin/env python

"""
A CLI tool for vetting sequences and sequence alignments.

The program will calculate pairwise sequence distances and rank sequences based
on their mean and maximum distances, which can be useful for identifying
potentially erroneous sequences. Also, any sequences involved in comparisons
that have a smaller distance when reverse and complemented will be reported.
"""

import os
import sys
import argparse

from seqsift.utils import FILE_FORMATS, VALID_DATA_TYPES, argparse_utils

_program_info = {
    'name': os.path.basename(__file__),
    'author': 'Jamie Oaks',
    'version': 'Version 0.1.0',
    'description': __doc__,
    'copyright': 'Copyright (C) 2013 Jamie Oaks',
    'license': 'GNU GPL version 3 or later',}

def main_cli():
    description = '{name} {version}'.format(**_program_info)
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument('input_file', metavar='INPUT-SEQ-FILE',
            type = argparse_utils.arg_is_file,
            help = ('Input sequence file to be vetted.'))
    parser.add_argument('-n', '--num-samples',
            type = int,
            default = 0,
            help = ('The number of randomly sampled sequences to which each '
                    'sequence will be compared. If less than 1 (the defualt is '
                    '0), all pairwise comparisons will be performed. For very '
                    'large numbers of sequences, performing all pairwise '
                    'comparisons will take a long time. This option will speed '
                    'things up as long as the number specified is less than '
                    'about half of the number of input sequences. If the '
                    'number you are considering is close to half of the number '
                    'sequences, you should probably specify zero and do all '
                    'combinations. You should not specify a number greater '
                    'half the number of sequences, because it will take longer '
                    'and be less thorough than the default.'))
    parser.add_argument('-a', '--aligned',
            action = 'store_true',
            help = ('Treat input sequences as aligned. I.e., do not perform '
                    'pairwise alignment before calculating distances between '
                    'sequences.'))
    parser.add_argument('--msa',
            action = 'store_true',
            help = ('Perform a full multiple sequence alignemnt prior to '
                    'comparing sequences. The default is to align each '
                    'pair of sequences being compared. This option is '
                    'overruled by the `-a`/`--aligned` option, which prevents '
                    'any aligning of the sequences. If this option is used '
                    'the resulting alignment is written to file.'))
    parser.add_argument('--aligner',
            type = argparse_utils.arg_is_executable,
            help = ('Path to alignment program executable to use for aligning. '
                    'The default is to look for mafft and then muscle in PATH, '
                    'and if neither are found use the (slow) built-in '
                    'function. The aligner is not used if the `-a`/`--aligned` '
                    'option is specified.'))
    parser.add_argument('--format',
            dest = 'input_format',
            type = str,
            choices = FILE_FORMATS.supported_formats,
            help = ('The format of the input sequence file. Valid options '
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
    parser.add_argument('-o', '--output-dir',
            type = argparse_utils.arg_is_dir,
            help = ('The directory in which all output files will be written. '
                    'The default is to use the directory of the input file.'))
    parser.add_argument('--seed',
            action = 'store',
            type = int,
            help = ('Random number seed to use for the analysis. This option '
                    'is only revelant if a number greater than 0 is specified '
                    'for the `-n/--num-samples` option.'))
    parser.add_argument('--quiet',
            action = 'store_true',
            help = 'Run with verbose messaging.')
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
    log = get_logger()

    ##########################################################################
    ## package imports

    from seqsift.utils import GLOBAL_RNG, dataio, functions
    from seqsift.seqops import seqstats
    from seqsift.utils.fileio import OpenFile

    ##########################################################################
    ## handle args

    ## set seed if randomly sampling sequences
    if args.num_samples > 0:
        if not args.seed:
            args.seed = random.randint(1, 999999999)
        GLOBAL_RNG.seed(args.seed)
        log.warning('Seed: {0}'.format(args.seed))

    ## get input file format
    if not args.input_format:
        args.input_format = FILE_FORMATS.get_format_from_file_object(
                args.input_file)
    if not args.input_format:
        log.error("Could not determine input format.\n"
                   "You must either provide the input format\n"
                   "using the '--from' option or have a recognizable\n"
                   "file extension on the input file name.\n"
                   "Here are the supported file extensions:\n{0}".format(
                        str(FILE_FORMATS)))
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)
    
    aligner_tools = ['mafft', 'muscle']
    if args.aligner:
        aligner_tools = [args.aligner]

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.input_file)

    full_alignment_out_path = os.path.join(
                args.output_dir, 'seqvet-msa.txt')

    ##########################################################################
    ## heavy lifting

    seqs = dataio.get_seq_iter(args.input_file,
            format = args.input_format,
            data_type = args.data_type)

    log.info('Calculating pairwise distances...')
    distances, rev_comp_errors = seqstats.summarize_distances(seqs,
            sample_size = args.num_samples,
            per_site = False,
            aligned = args.aligned,
            ignore_gaps = True,
            do_full_alignment = args.msa,
            full_alignment_out_path = full_alignment_out_path,
            aligner_tools = aligner_tools)
    log.info('Done!')

    log.info('Writing mean distances to file...')
    distances = sorted([(k, v) for k, v in distances.iteritems()],
            key = lambda x: x[1].mean,
            reverse = True)
    mean_path = functions.get_new_path(os.path.join(args.output_dir, 
            'seqvet-mean-distances.txt'))
    with OpenFile(mean_path, 'w') as out:
        out.write('seq_id\tmean_distance\n')
        for (seq_id, dist) in distances:
            out.write('{0}\t{1}\n'.format(seq_id, dist.mean))

    log.info('Writing max distances to file...')
    distances = sorted(distances,
            key = lambda x: x[1].maximum,
            reverse = True)
    max_path = functions.get_new_path(os.path.join(args.output_dir, 
            'seqvet-max-distances.txt'))
    with OpenFile(max_path, 'w') as out:
        out.write('seq_id\tmax_distance\n')
        for (seq_id, dist) in distances:
            out.write('{0}\t{1}\n'.format(seq_id, dist.maximum))

    if rev_comp_errors:
        rce = [sorted([s1, s2]) for (s1, s2, d, drc) in rev_comp_errors]
        rce = sorted(set([(s1, s2) for (s1, s2) in rce]))
        log.info('Writing potential reverse-complement errors to file...')
        path = functions.get_new_path(os.path.join(args.output_dir, 
                'seqvet-reverse-complement-warnings.txt'))
        with OpenFile(path, 'w') as out:
            out.write('seq1\tseq2\n')
            for (seq1, seq2) in rce:
                out.write('{0}\t{1}\n'.format(seq1, seq2))

if __name__ == '__main__':
    main_cli()

