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
import random
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

    comparison_args = parser.add_argument_group('Comparison Options',
            'Options to control the number and nature of sequence comparisons')
    comparison_args.add_argument('-n', '--num-samples',
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
                    'combinations. You should not specify a number greater than '
                    'half the number of sequences, because it will take longer '
                    'and be less thorough than the default.'))
    comparison_args.add_argument('--seed',
            action = 'store',
            type = int,
            help = ('Random number seed to use for the analysis. This option '
                    'is only revelant if a number greater than 0 is specified '
                    'for the `-n/--num-samples` option.'))
    comparison_args.add_argument('--compare-translated',
            action = 'store_true',
            help = ('Compare amino acid sequences encoded by the longest '
                'reading frame found in each sequence. To use this option, '
                '`data-type` must be dna or rna. See "Translation Options" '
                'for controlling how the longest reading frame of each '
                'sequence is determined and translated.'))
    comparison_args.add_argument('--check-ids',
            action = 'store_true',
            help = ('Check sequence IDs for duplicates.'))
    comparison_args.add_argument('--summarize-reading-frame-lengths',
            action = 'store_true',
            help = ('Report the length of the longest reading frame of '
                    'each sequence. See "Translation Options" for controlling '
                    'how reading frames are determined.'))
    comparison_args.add_argument('-g', '--count-gaps',
            action = 'store_true',
            help = ('Count gaps when calculating pairwise sequence distances. '
                    'The default is to calculate (number of differences '
                    'ignoring gaps / number of aligned sites ignoring sites '
                    'with gaps) for each pairwise comparison. When this option '
                    'is used, the distance is (number of differences including '
                    'gap differences / total number of aligned sites).'))

    alignment_args = parser.add_argument_group('Alignment Options',
            ('These options control if/how sequences are to be aligned prior '
             'to calculating distances.'))
    alignment_args.add_argument('-a', '--aligned',
            action = 'store_true',
            help = ('Treat input sequences as aligned. I.e., do not perform '
                    'pairwise alignment before calculating distances between '
                    'sequences (except when calculating distances for reverse '
                    'and complemented sequences).'))
    alignment_args.add_argument('--aligner',
            type = argparse_utils.arg_is_executable,
            help = ('Path to alignment program executable to use for pairwise'
                    'alignments of sequences. '
                    'The default is to look for muscle and then mafft in PATH, '
                    'and if neither are found use the (slow) built-in '
                    'function. Even if the `-a`/`--aligned` option is '
                    'specified, the aligner will still be used for pairwise '
                    'alignments when calculating distances of reverse and '
                    'complemented sequences.'))
    alignment_args.add_argument('--msa',
            action = 'store_true',
            help = ('Perform a full multiple sequence alignemnt prior to '
                    'comparing sequences. The default is to align each '
                    'pair of sequences being compared. This option is '
                    'overruled by the `-a`/`--aligned` option. '
                    'If this option is used '
                    'the resulting alignment is written to file.'))
    alignment_args.add_argument('--msa-aligner',
            type = argparse_utils.arg_is_executable,
            help = ('Path to alignment program executable to use for full '
                    'multiple sequence alignment. '
                    'The default is to look for mafft and then muscle in PATH, '
                    'and if neither are found the program will exit with an '
                    'error message. If you do not have mafft or muscle '
                    'you cannot use this option. '
                    'This option is only used if the `-a`/`--aligned` option '
                    'is not specified, and the `--msa` option is specified.'))

    translation_args = parser.add_argument_group('Translation Options',
            ('These options control translation from nucleotide to amino acid '
             'sequences.'))
    translation_args.add_argument('--table',
            type = int,
            choices = list(range(1, 7)) + list(range(9, 17)) + list(range(21, 26)),
            default = 1,
            help = ('The translation table to use for any options associated '
                    'with translating nucleotide sequences to amino acids. '
                    'Option should be the integer that corresponds to the '
                    'desired translation table according to NCBI '
                    '(http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). '
                    'The default is 1 (the "standard" code).'))
    translation_args.add_argument('--allow-partial',
            action = 'store_true',
            default = False,
            help = ('Allow partial reading frames at the beginning (no start '
                'codon) and end (no stop codon) of sequences.'))
    translation_args.add_argument('--read-after-stop',
            action = 'store_true',
            default = False,
            help = ('A new reading frame begins immediately after a stop codon. '
                    'The default is to start reading frame at next start codon '
                    'after a stop codon. This option might be useful for exons.'))

    data_args = parser.add_argument_group('Data Options',
            ('Options specifying the input data type and format'))
    data_args.add_argument('-d', '--data-type',
            type = str,
            choices = VALID_DATA_TYPES,
            default='dna',
            help = ('The type of sequence data. The default is dna. Valid '
                    'options include: {0}.'.format(', '.join(
                            VALID_DATA_TYPES))))
    data_args.add_argument('--format',
            dest = 'input_format',
            type = str,
            choices = FILE_FORMATS.supported_formats,
            help = ('The format of the input sequence file. Valid options '
                    'include: {0}. By default, the format is guessed based on '
                    'the extension of the first input file. However, if '
                    'provided, this option will always take precedence over '
                    'the file extension.'.format(
                          ', '.join(FILE_FORMATS.supported_formats))))

    output_args = parser.add_argument_group('Output Options',
            'Options for controlling output of program')
    output_args.add_argument('-o', '--output-dir',
            type = argparse_utils.arg_is_dir,
            help = ('The directory in which all output files will be written. '
                    'The default is to use the directory of the input file.'))

    messaging_args = parser.add_argument_group('Messaging Options',
            ('These options control verbosity of messaging.'))
    messaging_args.add_argument('--log-frequency',
            type = argparse_utils.arg_is_nonnegative_int,
            default = 1000,
            help = ('The frequency at which to log progress. Default is to log '
                    'every 1000 sequence comparisons.'))
    messaging_args.add_argument('--quiet',
            action = 'store_true',
            help = 'Run without verbose messaging.')
    messaging_args.add_argument('--debug',
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

    from seqsift.utils import GLOBAL_RNG, dataio, functions, alphabets
    from seqsift.seqops import seqsum, seqmod, seqstats
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
    
    aligner_tools = ['muscle', 'mafft']
    if args.aligner:
        aligner_tools = [args.aligner]
    full_aligner_tools = ['mafft', 'muscle']
    if args.msa_aligner:
        full_aligner_tools = [args.msa_aligner]

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.input_file)

    full_alignment_out_path = os.path.join(
                args.output_dir, 'seqvet-msa.txt')
    alphabet = alphabets.DnaAlphabet()
    if args.data_type in ['aa', 'protein']:
        alphabet = alphabets.ProteinAlphabet()

    if (args.summarize_reading_frame_lengths and 
            (not args.data_type in ['dna', 'rna'])):
        log.error("`--summarize-reading-frame-lengths` is only compatible "
                   "with DNA or RNA.")
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)

    if (args.compare_translated and (not args.data_type in ['dna', 'rna'])):
        log.error("`-compare-translated` is only compatible with DNA or RNA.")
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)

    ##########################################################################
    ## heavy lifting

    seqs = dataio.get_seq_iter([args.input_file],
            format = args.input_format,
            data_type = args.data_type)

    if args.summarize_reading_frame_lengths:
        log.info('Summarizing longest reading frame lengths...')
        if not isinstance(seqs, dataio.BufferedIter):
            seqs = dataio.BufferedIter(seqs)
        lengths = seqsum.summarize_longest_read_lengths(seqs,
                table = args.table,
                allow_partial = args.allow_partial,
                require_start_after_stop = (not args.read_after_stop))
        length_path = os.path.join(
                args.output_dir, 'seqvet-reading-frame-lengths.txt')
        log.info('Writing longest reading frame lengths to file...')
        with OpenFile(length_path, 'w') as out:
            out.write('seq_id\tlrf\trev_comp_lrf\n')
            for (l, rc_l, seq_id) in lengths:
                out.write('{0}\t{1}\t{2}\n'.format(seq_id, l, rc_l))

    if args.compare_translated:
        log.info('Translating longest reading frames for distance '
                 'calculations...')
        seqs = seqmod.translate_longest_reading_frames(seqs,
                table = args.table,
                allow_partial = args.allow_partial,
                require_start_after_stop = (not args.read_after_stop))
        alphabet = alphabets.ProteinAlphabet()

    if args.check_ids:
        log.info('Checking sequence IDs...')
        if not isinstance(seqs, dataio.BufferedIter):
            seqs = dataio.BufferedIter(seqs)
            dups = seqstats.get_duplicate_ids(seqs)
            if len(dups) > 0:
                dup_path = functions.get_new_path(os.path.join(args.output_dir,
                        'seqvet-duplicate-ids.txt'))
                log.warning('Duplicate IDs found! Writing them to '
                        '{0}'.format(dup_path))
                with OpenFile(dup_path, 'w') as out:
                    for dup in dups:
                        out.write('{0}\n'.format(dup))
            else:
                log.info('No duplicate sequence IDs were found.')

    log.info('Calculating pairwise distances...')
    distances, rev_comp_errors = seqsum.summarize_distances(seqs,
            sample_size = args.num_samples,
            per_site = True,
            aligned = args.aligned,
            ignore_gaps = (not args.count_gaps),
            alphabet = alphabet,
            do_full_alignment = args.msa,
            full_alignment_out_path = full_alignment_out_path,
            aligner_tools = aligner_tools,
            full_aligner_tools = full_aligner_tools,
            log_frequency = args.log_frequency)
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

