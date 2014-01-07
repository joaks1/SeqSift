#! /usr/bin/env python

import os
import shutil
from cStringIO import StringIO

from Bio.Align.Applications import MafftCommandline, MuscleCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqsift.utils.alphabets import DnaAlphabet
from seqsift.utils.fileio import TemporaryFilePath, OpenFile, expand_path
from seqsift.utils import functions, dataio, errors
from seqsift.seqops import sequtils
from seqsift.utils.messaging import get_logger

_LOG = get_logger(__name__)

DNA_AMBIGUITY_CODES = DnaAlphabet().residue_ambiguity_codes

def align(seq_iter, tools = ['mafft', 'muscle'], out_path = None):
    """
    Returns BufferedIter of aligned sequences in `seq_iter`. The
    `tools` argument should be a prioritized list of options for the external
    alignment program to use. For example, if `['mafft', 'muscle']` is
    specified (the default), mafft will be used if the executable is found in
    PATH. If mafft cannot be found, it will try to use muscle. If none of the
    listed programs can be found (or if the argument is an empty list or None),
    the (much slower) built-in `global_align` function is used.
    """
    aligner = get_aligner(tools = tools, out_path = out_path)
    if not aligner:
        raise errors.ExternalToolNotFoundError('No alignment tools found '
                'to perform multiple sequence alignment')
    return aligner.align(seq_iter)

def align_pair(seq_record1, seq_record2, tools = ['mafft', 'muscle']):
    """
    Returns the aligned copies of `SeqRecord` `seq_record1` and `seq_record2`.
    The `tools` argument should be a prioritized list of options for the
    external alignment program to use. For example, if `['mafft', 'muscle']` is
    specified (the default), mafft will be used if the executable is found in
    PATH. If mafft cannot be found, it will try to use muscle. If none of the
    listed programs can be found (or if the argument is an empty list or None),
    the (much slower) built-in `global_align` function is used.
    """
    aligner = get_aligner(tools = tools)
    if not aligner:
        if tools:
            _LOG.warning('WARNING: external alignment tools not found; '
                    'using (slow) built-in alignment function.')
        seq1, seq2 = global_align(seq_record1, seq_record2)
        s1 = sequtils.copy_seq_metadata(seq_record1, seq1)
        s2 = sequtils.copy_seq_metadata(seq_record2, seq2)
        return s1, s2
    seqs = list(aligner.align([seq_record1, seq_record2]))
    print seqs[0].seq
    print seqs[1].seq
    assert len(seqs) == 2
    sequences = dict(zip([s.id for s in seqs], seqs))
    return sequences[seq_record1.id], sequences[seq_record2.id]

def get_aligner(tools = ['mafft', 'muscle'], out_path = None):
    tool_dict = {}
    if tools:
        tool_dict = dict(zip([os.path.basename(t).lower() for t in tools],
                tools))
    for tool, exe in tool_dict.iteritems():
        if tool == 'mafft':
            try:
                return MafftAligner(exe = exe, out_path = out_path)
            except errors.ExternalToolNotFoundError as e:
                _LOG.warning('WARNING: mafft is not available')
        elif tool == 'muscle':
            try:
                return MuscleAligner(exe = exe, out_path = out_path)
            except errors.ExternalToolNotFoundError as e:
                _LOG.warning('WARNING: muscle is not available')
        else:
            _LOG.warning('WARNING: {0} is not a supported alignment '
                    'tool'.format(tool))
    return None

def global_align(
        seq1,
        seq2,
        similarity_matrix = {
            ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
        },
        gap_cost = -5,
        gap_char = '-'):
    """
    Returns the global pairwise alignment of `seq1` and `seq2` using the
    Needleman-Wunsch algorithm.
    """
    seq1 = [x for x in seq1 if x != gap_char]
    seq2 = [x for x in seq2 if x != gap_char]
    fmatrix = calculate_F_matrix(
            seq1=seq1, seq2=seq2,
            similarity_matrix=similarity_matrix,
            gap_cost=gap_cost)
    s1, s2 = trace_max_score(fmatrix=fmatrix, seq1=seq1, seq2=seq2,
            similarity_matrix=similarity_matrix,
            gap_cost=gap_cost,
            gap_char=gap_char)
    return s1, s2

def calculate_F_matrix(
        seq1,
        seq2,
        similarity_matrix = {
            ('A', 'A'): 10, ('A', 'G'): -1, ('A', 'C'): -3, ('A', 'T'): -4,
            ('G', 'A'): -1, ('G', 'G'):  7, ('G', 'C'): -5, ('G', 'T'): -3,
            ('C', 'A'): -3, ('C', 'G'): -5, ('C', 'C'):  9, ('C', 'T'):  0,
            ('T', 'A'): -4, ('T', 'G'): -3, ('T', 'C'):  0, ('T', 'T'):  8,
        },
        gap_cost = -5):
    """
    Computes and returns the F matrix of the Needleman-Wunsch algorithm for
    `seq1` and `seq2`.
    """
    nrows = len(seq1) + 1
    ncols = len(seq2) + 1
    fmatrix = [[None for j in range(ncols)] for i in range(nrows)]
    for i in range(nrows):
        fmatrix[i][0] = i * gap_cost
    for j in range(ncols):
        fmatrix[0][j] = j * gap_cost
    for i in range(1, nrows):
        for j in range(1, ncols):
            match = fmatrix[i - 1][j - 1] + \
                        match_score(
                                seq1[i - 1].upper(),
                                seq2[j - 1].upper(),
                                similarity_matrix)
            deletion = fmatrix[i - 1][j] + gap_cost
            insertion = fmatrix[i][j - 1] + gap_cost
            fmatrix[i][j] = max(match, deletion, insertion)
    return fmatrix

def match_score(base1, base2, similarity_matrix):
    if (base1, base2) in similarity_matrix.keys():
        return similarity_matrix[base1, base2]
    elif (base1 in DNA_AMBIGUITY_CODES.keys()) or \
            (base2 in DNA_AMBIGUITY_CODES.keys()):
        possible_bases1 = DNA_AMBIGUITY_CODES.get(base1, base1)
        possible_bases2 = DNA_AMBIGUITY_CODES.get(base2, base2)
        total_score = 0
        comparisons = 0
        for b1 in possible_bases1:
            for b2 in possible_bases2:
                total_score += similarity_matrix[b1, b2]
                comparisons += 1
        assert comparisons == len(possible_bases1) * len(possible_bases2)
        return float(total_score)/comparisons
    else:
        raise ValueError(
                "{0!r} or {1!r} was not found ".format(base1, base2) + \
                "in the similarity matrix or ambiguity codes.")

def trace_max_score(
        fmatrix,
        seq1,
        seq2,
        similarity_matrix,
        gap_cost,
        gap_char = '-'):
    """
    Returns the global pairwise alignment of `seq1` and `seq2` using the
    Needleman-Wunsch algorithm, given the `fmatrix`.
    """
    s1 = ''
    s2 = ''
    i = len(seq1)
    j = len(seq2)
    while i > 0 and j > 0:
        score = fmatrix[i][j]
        diagonal_score = fmatrix[i - 1][j - 1]
        up_score = fmatrix[i - 1][j]
        left_score = fmatrix[i][j - 1]
        if score == diagonal_score + match_score(
                seq1[i - 1].upper(),
                seq2[j - 1].upper(),
                similarity_matrix):
            s1 = seq1[i - 1] + s1
            s2 = seq2[j - 1] + s2
            i -= 1
            j -= 1
        elif score == up_score + gap_cost:
            s1 = seq1[i - 1] + s1
            s2 = gap_char + s2
            i -= 1
        elif score == left_score + gap_cost:
            s1 = gap_char + s1
            s2 = seq2[j - 1] + s2
            j -= 1
        else:
            raise Exception("Error in tracing F Matrix.")
    while i > 0:
        s1 = seq1[i - 1] + s1
        s2 = gap_char + s2
        i -= 1
    while j > 0:
        s1 = gap_char + s1
        s2 = seq2[j - 1] + s2
        j -= 1
    return s1, s2

class MafftAligner(object):
    count = 0
    def __init__(self, exe = 'mafft', out_path = None, **kwargs):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.exe = functions.which(exe)
        if not self.exe:
            raise errors.ExternalToolNotFoundError(
                    'Cannot find mafft executable')
        _LOG.info('{0}: Using exe {1!r}'.format(self.name, self.exe))
        self.kwargs = kwargs
        if self.kwargs.has_key('input'):
            raise ValueError('MafftAligner does not accept the keyword '
                    'argument `input`.')
        if not self.kwargs:
            self.kwargs = {'auto': True}
        self.out_path = out_path
        self.cmd = None
    
    def align(self, seq_iter):
        with TemporaryFilePath() as tmp_path:
            with OpenFile(tmp_path, 'w') as tmp:
                for seq in seq_iter:
                    tmp.write(seq.format('fasta'))
            self.kwargs['input'] = tmp_path
            mafft_command = MafftCommandline(self.exe, **self.kwargs)
            self.cmd = str(mafft_command)
            _LOG.debug('{0}: Executing command {1!r}'.format(self.name,
                    self.cmd))
            stdout, stderr = mafft_command()
            if self.out_path:
                self.out_path = functions.get_new_path(self.out_path)
                with OpenFile(self.out_path, 'w') as out:
                    out.write(stdout)
        return dataio.get_buffered_seq_iter(StringIO(stdout), format='fasta')

class MuscleAligner(object):
    count = 0
    def __init__(self, exe = 'muscle', out_path = None, **kwargs):
        self.__class__.count += 1
        self.name = self.__class__.__name__ + '-' + str(self.count)
        self.exe = functions.which(exe)
        if not self.exe:
            raise errors.ExternalToolNotFoundError(
                    'Cannot find muscle executable')
        _LOG.info('{0}: Using exe {1!r}'.format(self.name, self.exe))
        self.kwargs = kwargs
        if self.kwargs.has_key('input') or self.kwargs.has_key('out'):
            raise ValueError('MuscleAligner does not accept keyword '
                    'arguments `input`/`out`.')
        if out_path:
            out_path = expand_path(out_path)
        self.out_path = out_path
        self.cmd = None
    
    def align(self, seq_iter):
        with TemporaryFilePath() as in_path:
            with OpenFile(in_path, 'w') as tmp:
                for seq in seq_iter:
                    tmp.write(seq.format('fasta'))
            self.kwargs['input'] = in_path
            with TemporaryFilePath() as tmp_out_path:
                self.kwargs['out'] = tmp_out_path
                muscle_command = MuscleCommandline(self.exe, **self.kwargs)
                self.cmd = str(muscle_command)
                _LOG.debug('{0}: Executing command {1!r}'.format(self.name,
                        self.cmd))
                stdout, stderr = muscle_command()
                results = dataio.get_buffered_seq_iter(tmp_out_path,
                        format='fasta')
                if self.out_path:
                    self.out_path = functions.get_new_path(self.out_path)
                    shutil.mv(tmp_out_path, self.out_path)
        return results

