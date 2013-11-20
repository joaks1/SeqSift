#! /usr/bin/env python

from seqsift.utils.alphabets import DnaAlphabet

DNA_AMBIGUITY_CODES = DnaAlphabet().residue_ambiguity_codes

def global_align(
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
    Returns the global pairwise alignment of `seq1` and `seq2` using the
    Needleman-Wunsch algorithm.
    """
    fmatrix = calculate_F_matrix(
            seq1=seq1, seq2=seq2,
            similarity_matrix=similarity_matrix,
            gap_cost=gap_cost)
    s1, s2 = trace_max_score(fmatrix=fmatrix, seq1=seq1, seq2=seq2,
            similarity_matrix=similarity_matrix,
            gap_cost=gap_cost)
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
        gap_cost):
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
            s2 = '-' + s2
            i -= 1
        elif score == left_score + gap_cost:
            s1 = '-' + s1
            s2 = seq2[j - 1] + s2
            j -= 1
        else:
            raise Exception("Error in tracing F Matrix.")
    while i > 0:
        s1 = seq1[i - 1] + s1
        s2 = '-' + s2
        i -= 1
    while j > 0:
        s1 = '-' + s1
        s2 = seq2[j - 1] + s2
        j -= 1
    return s1, s2
