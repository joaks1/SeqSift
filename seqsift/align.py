#!/usr/bin/python

def global_align(
        seq1,
        seq2,
        similarity_matrix = {
        ('a', 'a'): 10, ('a', 'g'): -1, ('a', 'c'): -3, ('a', 't'): -4,
        ('g', 'a'): -1, ('g', 'g'):  7, ('g', 'c'): -5, ('g', 't'): -3,
        ('c', 'a'): -3, ('c', 'g'): -5, ('c', 'c'):  9, ('c', 't'):  0,
        ('t', 'a'): -4, ('t', 'g'): -3, ('t', 'c'):  0, ('t', 't'):  8,
        },
        gap_cost = -5):
    """
    Returns the global pairwise alignment of `seq1` and `seq2` using the
    Needleman-Wunsch algorithm.
    """
    fmatrix = computeFMatrix(
            seq1=seq1, seq2=seq2,
            similarity_matrix=similarity_matrix,
            gap_cost=gap_cost)
    s1, s2 = trace_max_score(fmatrix=fmatrix, seq1=seq1, seq2=seq2,
            similarity_matrix=similarity_matrix,
            gap_cost=gap_cost)
    return s1, s2

def computeFMatrix(
        seq1,
        seq2,
        similarity_matrix = {
        ('a', 'a'): 10, ('a', 'g'): -1, ('a', 'c'): -3, ('a', 't'): -4,
        ('g', 'a'): -1, ('g', 'g'):  7, ('g', 'c'): -5, ('g', 't'): -3,
        ('c', 'a'): -3, ('c', 'g'): -5, ('c', 'c'):  9, ('c', 't'):  0,
        ('t', 'a'): -4, ('t', 'g'): -3, ('t', 'c'):  0, ('t', 't'):  8,
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
                    similarity_matrix[seq1[i - 1].lower(), seq2[j - 1].lower()]
            deletion = fmatrix[i - 1][j] + gap_cost
            insertion = fmatrix[i][j - 1] + gap_cost
            fmatrix[i][j] = max(match, deletion, insertion)
    return fmatrix

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
        if score == diagonal_score + similarity_matrix[
                seq1[i - 1].lower(),
                seq2[j - 1].lower()]:
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
