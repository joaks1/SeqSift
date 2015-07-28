#! /usr/bin/env python

import os
import sys
import re
from cStringIO import StringIO

from seqsift.utils import dataio, stats, external_tools, tempfs, alphabets
from seqsift.utils.fileio import OpenFile
from seqsift.seqops import seqstats, seqfilter
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

number_pattern_string = r'[\d\.Ee\-\+]+'
paup_hky_string = (
        '^\s*Estimated\s+base\s+frequencies\s+=\s+'
        'A:(?P<a>{0})\s+C:(?P<c>{0})\s+'
        'G:(?P<g>{0})\s+T:(?P<t>{0})'
        '\n\s*Estimated\s+ti/tv\s+ratio\s+=\s+{0}\s*'
        '\(kappa\s+=\s+(?P<kappa>{0})\)\s*$'.format(number_pattern_string)
        )
paup_hky_pattern = re.compile(paup_hky_string, re.MULTILINE)

class PyMsBayesComparisons(object):
    count = 0
    alphabet = alphabets.DnaAlphabet()
    sample_table_header = ('# taxon\tlocus\tploidy_multiplier\t'
                'rate_multiplier\tnsamples1\tnsamples2\tkappa\t'
                'nsites\ta\tc\tg\tpath\n'
                'BEGIN SAMPLE_TBL\n')
    sample_table_footer = 'END SAMPLE_TBL\n'
    pi_header = 'taxon\tlocus\tpi1\tpi2\n'
    def __init__(self,
            comparisons,
            name = None,
            locus = None):
        self.__class__.count += 1
        if not name:
            name = self.__class__.__name__ + '-' + str(self.count)
        self.name = name
        if not locus:
            locus = 'locus{0}'.format(self.count)
        self.locus = locus
        self.comparisons = comparisons

    def add_sequence(self, seq_record, strict = False):
        ret = None
        for comp in self.comparisons:
            r = comp.add_sequence(seq_record)
            if r is not None:
                ret = r
        if strict and (ret is None):
            raise Exception('{0} is not in comparisons'.format(seq_record.id))

    def extend_sequences(self, seq_record_iter, strict = False):
        for s in seq_record_iter:
            self.add_sequence(s, strict = strict)

    def _get_smallest_number_of_sequences(self):
        smallest = None
        for comp in self.comparisons:
            s = comp.smallest_sample_size
            if smallest is None:
                smallest = s
            if s < smallest:
                smallest = s
        return smallest

    smallest_sample_size = property(_get_smallest_number_of_sequences)

    def _get_shortest_alignment(self):
        shortest = None
        for comp in self.comparisons:
            s = comp.alignment_length
            if shortest is None:
                shortest = s
            if s < shortest:
                shortest = s
        return shortest

    shortest_alignment = property(_get_shortest_alignment)
    
    def write_comparisons(self, fasta_dir,
            config_dir = None,
            estimate_hky_parameters = False):
        s = StringIO()
        pi_s = StringIO()
        for comp in self.comparisons:
            if estimate_hky_parameters and (not comp.estimated_hky_parameters):
                comp.estimate_hky_parameters()
            nsamples = comp.number_of_sequences
            taxon = comp.comparison_str
            path = os.path.join(fasta_dir, '{0}-{1}.fasta'.format(
                    taxon, self.locus))
            rel_path = path
            if config_dir:
                rel_path = os.path.relpath(path, config_dir)
            comp.write_sequences(path)
            s.write('{taxon}\t{locus}\t{ploidy_multiplier}\t'
                '{rate_multiplier}\t{nsamples1}\t{nsamples2}\t{kappa}\t'
                '{nsites}\t{a}\t{c}\t{g}\t{path}\n'.format(
                    taxon = taxon,
                    locus = self.locus,
                    ploidy_multiplier = comp.ploidy_multiplier,
                    rate_multiplier = comp.rate_multiplier,
                    nsamples1 = nsamples[0],
                    nsamples2 = nsamples[1],
                    kappa = comp.kappa,
                    nsites = comp.alignment_length,
                    a = comp.a,
                    c = comp.c,
                    g = comp.g,
                    path = rel_path))
            comp.estimate_pi()
            pi_s.write('{taxon}\t{locus}\t{pi1}\t{pi2}\n'.format(
                    taxon = taxon,
                    locus = self.locus,
                    pi1 = comp.pi[0],
                    pi2 = comp.pi[1]))
        return s.getvalue(), pi_s.getvalue()

    @classmethod
    def process_loci_file(cls,
            loci_file_obj,
            pop_id_maps,
            fasta_out_dir,
            config_out_dir = None,
            minimum_sample_size = 2,
            minimum_alignment_length = 50,
            max_ambiguities_per_seq = 0.2,
            require_shared_loci = False,
            estimate_hky_parameters = False):
        if not config_out_dir:
            config_out_dir = fasta_out_dir
        config_out_path = os.path.join(config_out_dir,
                    'config-sample-table.txt')
        pi_out_path = os.path.join(config_out_dir,
                    'pi-estimates.txt')
        config_stream = StringIO()
        pi_stream = StringIO()
        for i, locus in enumerate(dataio.LociFileIter(loci_file_obj)):
            comps = []
            for id_map in pop_id_maps:
                assert(len(id_map) == 2)
                pop1, pop2 = sorted(id_map.keys())
                c = Comparison(
                        population_1_ids = id_map[pop1],
                        population_2_ids = id_map[pop2],
                        population_1_name = pop1,
                        population_2_name = pop2,
                        )
                # remove rows with many ambiguities
                seqs = seqfilter.row_filter(locus,
                        character_list = (cls.alphabet.ambiguity_codes.keys() +
                                [cls.alphabet.gap]),
                        max_frequency = max_ambiguities_per_seq)
                # remove all columns with ambiguities
                seqs = seqfilter.column_filter(seqs,
                        character_list = (cls.alphabet.ambiguity_codes.keys() +
                                [cls.alphabet.gap]),
                        max_frequency = 0.00001)
                c.extend_sequences(seqs)
                if ((c.alignment_length >= minimum_alignment_length) and
                        (c.smallest_sample_size >= minimum_sample_size)):
                    comps.append(c)
            if not comps:
                continue
            if require_shared_loci and (len(comps) < len(pop_id_maps)):
                continue
            pymsbayes_comps = cls(comparisons = comps,
                    locus = 'locus{0}'.format(i))
            config_str, pi_str = pymsbayes_comps.write_comparisons(
                    fasta_dir = fasta_out_dir,
                    config_dir = config_out_dir,
                    estimate_hky_parameters = estimate_hky_parameters)
            config_stream.write(config_str)
            pi_stream.write(pi_str)
        with OpenFile(config_out_path, 'w') as out:

            out.write(cls.sample_table_header)
            out.write(config_stream.getvalue())
            out.write(cls.sample_table_footer)
        with OpenFile(pi_out_path, 'w') as out:
            out.write(cls.pi_header)
            out.write(pi_stream.getvalue())
        return i + 1

    @classmethod
    def process_loci_files_as_pairs(cls,
            loci_file_objects,
            id_component_delimiter,
            id_component_index,
            fasta_out_dir,
            config_out_dir = None,
            minimum_sample_size = 2,
            minimum_alignment_length = 50,
            max_ambiguities_per_seq = 0.2,
            estimate_hky_parameters = False):
        if not config_out_dir:
            config_out_dir = fasta_out_dir
        config_out_path = os.path.join(config_out_dir,
                    'config-sample-table.txt')
        pi_out_path = os.path.join(config_out_dir,
                    'pi-estimates.txt')
        config_stream = StringIO()
        pi_stream = StringIO()
        locus_count = 0
        for i, loci_file_obj in enumerate(loci_file_objects):
            comparison_prefix = None
            try:
                comparison_prefix = os.path.basename(
                        loci_file_obj).split('.')[0]
            except:
                try:
                    comparison_prefix = os.path.basename(
                            loci_file_obj.name).split('.')[0]
                except:
                    pass
                pass
            for j, locus in enumerate(dataio.LociFileIter(loci_file_obj)):
                locus_count += 1
                # remove rows with many ambiguities
                seqs = seqfilter.row_filter(locus,
                        character_list = (cls.alphabet.ambiguity_codes.keys() +
                                [cls.alphabet.gap]),
                        max_frequency = max_ambiguities_per_seq)
                # remove all columns with ambiguities
                seqs = seqfilter.column_filter(seqs,
                        character_list = (cls.alphabet.ambiguity_codes.keys() +
                                [cls.alphabet.gap]),
                        max_frequency = 0.00001)
                c = Comparison.get_comparison_from_seqs(
                        sequences = seqs,
                        id_component_delimiter = id_component_delimiter,
                        id_component_index = id_component_index,
                        comparison_prefix = comparison_prefix,
                        )
                if not ((c.alignment_length >= minimum_alignment_length) and
                        (c.smallest_sample_size >= minimum_sample_size)):
                    continue
                pymsbayes_comps = cls(comparisons = [c])
                config_str, pi_str = pymsbayes_comps.write_comparisons(
                        fasta_dir = fasta_out_dir,
                        config_dir = config_out_dir,
                        estimate_hky_parameters = estimate_hky_parameters)
                config_stream.write(config_str)
                pi_stream.write(pi_str)
        with OpenFile(config_out_path, 'w') as out:

            out.write(cls.sample_table_header)
            out.write(config_stream.getvalue())
            out.write(cls.sample_table_footer)
        with OpenFile(pi_out_path, 'w') as out:
            out.write(cls.pi_header)
            out.write(pi_stream.getvalue())
        return locus_count


class Comparison(object):
    count = 0
    def __init__(self,
            population_1_ids,
            population_2_ids,
            sequences = None,
            population_1_name = None,
            population_2_name = None,
            comparison_prefix = None,
            name = None):
        self.__class__.count += 1
        if not name:
            name = self.__class__.__name__ + '-' + str(self.count)
        self.name = name
        s1 = set(population_1_ids)
        s2 = set(population_2_ids)
        if len(s1) != len(population_1_ids):
            raise Exception('duplicate sequence ids for population 1')
        if len(s2) != len(population_2_ids):
            raise Exception('duplicate sequence ids for population 2')
        if len(s1.intersection(s2)) > 0:
            raise Exception('cannot have the same sequence id in both '
                    'populations')
        self.population_1_ids = population_1_ids
        self.population_2_ids = population_2_ids
        self.population_1 = Population(name = population_1_name)
        self.population_2 = Population(name = population_2_name)
        self.comparison_str = '{0}-{1}'.format(self.population_1.name,
                self.population_2.name)
        if comparison_prefix:
            self.comparison_str = comparison_prefix + '-' + self.comparison_str
        self.alignment_length = None
        self.pi = (None, None)
        self.ploidy_multiplier = 1.0
        self.rate_multiplier = 1.0
        self.kappa = 1.0
        self.a = 0.25
        self.c = 0.25
        self.g = 0.25
        self.estimated_hky_parameters = False
        if sequences:
            self.extend_sequences(sequences)

    @classmethod
    def get_comparison_from_seqs(cls,
            sequences,
            id_component_delimiter,
            id_component_index,
            comparison_prefix = None,
            comparison_name = None):
        seq_list = [s for s in sequences]
        pops = {}
        for seq in seq_list:
            pop_id = seq.id.split(id_component_delimiter)[id_component_index]
            if pop_id in pops:
                pops[pop_id].append(seq.id)
            else:
                pops[pop_id] = [seq.id]
        if len(pops) != 2:
            raise Exception('Problem parsing populations from sequences:\n'
                    '{0} populations found ({1})'.format(len(pops),
                            ', '.join(pops.keys())))
        pop_1_id, pop_2_id = pops.keys()
        return cls(population_1_ids = pops[pop_1_id],
                population_2_ids = pops[pop_2_id],
                sequences = seq_list,
                population_1_name = pop_1_id,
                population_2_name = pop_2_id,
                comparison_prefix = comparison_prefix,
                name = comparison_name)

    def add_sequence(self, seq_record, strict = False):
        if seq_record.id in self.population_1_ids:
            self.population_1.add_sequence(seq_record)
            if self.alignment_length is None:
                self.alignment_length = len(seq_record.seq)
            if not (seq_record.seq) != self.alignment_length:
                raise Exception('sequence {0} is not length {1}: sequences must be '
                        'aligned'.format(seq_record.id, self.alignment_length))
            return seq_record.id
        elif seq_record.id in self.population_2_ids:
            self.population_2.add_sequence(seq_record)
            if self.alignment_length is None:
                self.alignment_length = len(seq_record.seq)
            if not (seq_record.seq) != self.alignment_length:
                raise Exception('sequence {0} is not length {1}: sequences must be '
                        'aligned'.format(seq_record.id, self.alignment_length))
        return seq_record.id
        if strict:
            raise Exception('{0} not in population ids'.format(seq_record.id))
        else:
            return None

    def extend_sequences(self, seq_record_iter, strict = False):
        for s in seq_record_iter:
            self.add_sequence(s, strict = strict)

    def write_sequences(self, dest = None):
        keys1 = sorted(self.population_1.sequences.keys())
        keys2 = sorted(self.population_2.sequences.keys())
        dataio.write_seqs(seqs = (
                [self.population_1.sequences[k] for k in keys1] + 
                [self.population_2.sequences[k] for k in keys2]),
                dest = dest,
                format = 'fasta')

    def _get_number_of_sequences(self):
        return (self.population_1.number_of_sequences,
                self.population_2.number_of_sequences)

    number_of_sequences = property(_get_number_of_sequences)

    def _get_smallest_number_of_sequences(self):
        return min(self.number_of_sequences)

    smallest_sample_size = property(_get_smallest_number_of_sequences)

    def estimate_pi(self):
        self.pi = (seqstats.average_number_of_pairwise_differences(
                seq_iter = self.population_1.sequences.values(),
                aligned = True,
                per_site = True),
             seqstats.average_number_of_pairwise_differences(
                seq_iter = self.population_2.sequences.values(),
                aligned = True,
                per_site = True))

    def estimate_hky_parameters(self):
        temp_fs = tempfs.TempFileSystem(parent = '/tmp')
        tmp_fasta_path = temp_fs.get_file_path()
        tmp_nex_data_path = temp_fs.get_file_path()
        tmp_nex_exe_path = temp_fs.get_file_path()
        tmp_score_path = temp_fs.get_file_path()
        self.write_sequences(dest = tmp_fasta_path)
        nseqs = dataio.convert_format(in_file = tmp_fasta_path,
                out_file = tmp_nex_data_path,
                in_format = 'fasta',
                out_format = 'nexus')
        write_paup_hky_file(tmp_nex_exe_path, tmp_nex_data_path,
                tmp_score_path)
        stdout_path = temp_fs.get_file_path()
        pw = external_tools.PaupWorker(nex_path = tmp_nex_exe_path,
                stdout_path = stdout_path)
        pw.start()
        parameters = parse_paup_log_file(stdout_path)
        self.kappa = parameters['kappa']
        if (self.kappa < 1.0) or (self.kappa > 1000.0):
            self.kappa = 1.0
        self.a = parameters['a']
        self.c = parameters['c']
        self.g = parametets['g']
        self.estimated_hky_parameters = True

class Population(object):
    count = 0
    def __init__(self,
            sequences = None,
            name = None):
        self.__class__.count += 1
        if not name:
            name = self.__class__.__name__ + '-' + str(self.count)
        self.name = name
        self.sequences = {}
        self.alignment_length = None
        if sequences:
            self.extend_sequences(sequences)

    def add_sequence(self, seq_record):
        if seq_record.id in self.sequences:
            raise Exception('{0} is already in {1}'.format(seq_record.id,
                    self.name))
        if self.alignment_length is None:
            self.alignment_length = len(seq_record.seq)
        if not (seq_record.seq) != self.alignment_length:
            raise Exception('sequence {0} is not length {1}: sequences must be '
                    'aligned'.format(seq_record.id, self.alignment_length))
        self.sequences[seq_record.id] = seq_record

    def extend_sequences(self, seq_record_iter, strict = False):
        for s in seq_record_iter:
            self.add_sequence(s, strict = strict)

    def _get_number_of_sequences(self):
        return len(self.sequences)

    number_of_sequences = property(_get_number_of_sequences)



def write_paup_hky_file(path, nexus_data_path, score_file_path):
    with open(path, 'w') as out:
        out.write("#NEXUS\n\n"
                "BEGIN PAUP;\n\t"
                "set warnreset=no; set increase=auto; set warnroot=no;\n\t"
                "execute {0};\n\t"
                "DSet distance=LOGDET objective=ME base=equal rates=equal pinv=0\n\tsubst=all negbrlen=setzero;\n\t"
                "NJ showtree=no breakties=random;\n"
                "END;\n\n"
                "BEGIN PAUP;\n\t"
                "Default lscores longfmt=yes;\n\t"
                "Set criterion=like;\n\t"
                "[!** Calculating HKY **]\n\t"
                "lscores 1/ nst=2 base=est tratio=est rates=equal pinv=0\n\t"
                "scorefile={1} append;\n"
                "END;\n".format(nexus_data_path, score_file_path)
                )

def parse_paup_log_file(path):
    with open(path, 'r') as in_stream:
        log = in_stream.read()
        match_iter = paup_hky_pattern.finditer(log)
        i = -1
        for i, m in enumerate(match_iter):
            pass
        if i < 0:
            return {'kappa': 1.0,
                    'a': 0.25,
                    'c': 0.25,
                    'g': 0.25,
                    't': 0.25}
        elif i > 0:
            raise Exception('found multiple sets of HKY parameters in '
                    '{0!r}'.format(path))
    return m.groupdict()


