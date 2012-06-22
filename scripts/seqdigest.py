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

from seqsift.utils import mkdr
from seqsift.digest import Fragment, RecognitionSeq, DigestSummary
from seqsift.utils.dataio import fetch_gb_seqs, get_seq_iter
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

_LOG = get_logger(__name__, 'INFO')

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

def digest_seq(recognition_seq, seq_record, extra_length=0,
               include_overhang=True):
    _LOG.info("Digesting seq {0} with recognition seq {1}...".format(
            seq_record.id, str(recognition_seq.seq)))
    return DigestSummary(recognition_seq=recognition_seq,
                         seq_record=seq_record,
                         extra_length=extra_length,
                         include_overhang=include_overhang)
               
def main():
    description = '{name} {version}'.format(**_program_info)
    usage = ("\n  %prog [options] -s <RECOGNITION_SEQUENCE> [<GENBANK_FILE1> "
             "<GENBANK_FILE2> ...]")
    parser = OptionParser(usage=usage, description=description,
                          version=_program_info['version'],
                          add_help_option=True)
    parser.add_option("-s", "--recognition_seq", dest="recognition_seq",
        type="string",
        help="Recognition sequence, 5' to 3', of restriction enzyme.")
    parser.add_option("-c", "--cut_site", dest="cut_site",
        type="int",
        help=("One-based index of the last base before the cut site in the "
              "recognition sequence. E.g., NotI: 5'---GC \ GGCCGC---3' "
              "has a cut site of 2, and would be passed to this program "
              "with '-s GCGGCCGC -c 2'."))
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
    parser.add_option("-x", "--extra_length", dest="extra_length", type="int",
        default=0,
        help=("Extra length (in bases) to add to each fragment. For example, "
              "you can include the length of oligos ligated to each fragment. "
              "This length is only added once, so if you want to simulate the "
              "ligation of oligos to both ends, of each fragment, provide the "
              "TOTAL length of the oligos ligated to each fragment. Do not "
              "add length for the overhang left by the restriction enzyme. "
              "That is handled by the program."))
    parser.add_option("--min_length", dest="min_length", type="int", default=0,
        help="Minimum fragment length to include in count.")
    parser.add_option("--max_length", dest="max_length", type="int",
        help="Maximum fragment length to include in count.")
    parser.add_option("-o", "--output_dir", dest="output_dir", type="string",
        help="Path to output directory. Default is './digests/'")
    (options, args) = parser.parse_args()
    
    if not options.recognition_seq or not options.cut_site:
        _LOG.error("You must provide a recognition sequence and cut site")
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)
    if not options.gi_numbers and len(args) < 1:
        _LOG.error("You must provide a sequence to digest, either via the "
                   "gi number option or sequence file arguments")
        sys.stderr.write(str(parser.print_help()))
        sys.exit(1)
    if options.format:
        format = EXTENSIONS[options.format.lower()]
    if not options.output_dir:
        out_dir = os.path.abspath(os.path.join(os.path.curdir, 'digests'))
        mkdr(out_dir)
    else:
        out_dir = os.path.expanduser(os.path.expandvars(options.output_dir))
    if not os.path.isdir(out_dir):
        _LOG.error("Output path {0} is not a directory".format(out_dir))
        sys.exit(1)
    if options.max_length and options.max_length < options.min_length:
        _LOG.error(
                "max_length ({0}) cannot be less than min_length ({1})".format(
                        options.max_length, options.min_length))
        sys.exit(1)
    if options.max_length:
        ml = str(options.max_length)
    else:
        ml = 'max'
        
    rseq = RecognitionSeq(str(options.recognition_seq), options.cut_site)
    gi_list = []
    if options.gi_numbers:
        gi_list = parse_gi_numbers(options.gi_numbers)
    
    digests = {}
    for gi in gi_list:
        _LOG.info("Downloading gi {0}...".format(gi))
        seq_iter = fetch_gb_seqs(str(gi), data_type='dna')
        for seq in seq_iter:
            digests[seq.id] = digest_seq(rseq, seq, options.extra_length)

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
            digests[seq.id] = digest_seq(rseq, seq, options.extra_length)
        f.close()

    _LOG.info('Finished digests!')
    _LOG.info('Summarizing results...')
    combined = {}
    total_count = 0
    filter_count = 0
    for id, digest in digests.iteritems():
        digest_total = 0
        digest_filter = 0
        out_file_path = os.path.join(out_dir, 
                ".".join([digest.molecule_name, 'txt']))
        out = open(out_file_path, 'w')
        out.write("{0}\t{1}\n".format('fragment_length', 'frequency'))
        for l in sorted(digest.length_distribution.iterkeys()):
            f = digest.length_distribution[l]
            out.write("{0}\t{1}\n".format(l, f))
            digest_total += f
            if options.max_length:
                if l <= options.max_length and l >= options.min_length:
                    digest_filter += f
            else:
                if l >= options.min_length:
                    digest_filter += f
            if l not in combined.keys():
                combined[l] = 0
            combined[l] += f
        out.close()
        _LOG.info('\nMolecule id: {0}\n'.format(digest.molecule_id) +
                  'Molecule name: {0}\n'.format(digest.molecule_name) +
                  'Molecule description: {0}\n'.format(
                        digest.molecule_description) +
                  '\ttotal fragments: {0}\n'.format(digest_total) +
                  '\tfragments of length {0}-{1}: {2}\n'.format(
                        options.min_length,
                        ml,
                        digest_filter) +
                  '\tfragment length distribution written to {0}\n'.format(
                        out_file_path))
        total_count += digest_total
        filter_count += digest_filter
    out_file_path = os.path.join(out_dir, 'combined.txt')
    out = open(out_file_path, 'w')
    out.write("{0}\t{1}\n".format('fragment_length', 'frequency'))
    for l in sorted(combined.iterkeys()):
        f = combined[l]
        out.write("{0}\t{1}\n".format(l, f))
    out.close()
    _LOG.info('\nSummary over ALL molecules:\n'
              '\ttotal fragments: {0}\n'.format(total_count) +
              '\tfragments of length {0}-{1}: {2}\n'.format(
                    options.min_length,
                    ml,
                    filter_count) +
              '\tfragment length distribution written to {0}\n'.format(
                    out_file_path))

if __name__ == '__main__':
    main()
