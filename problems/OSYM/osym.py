#!/usr/bin/env python
"""
Isolating Symbols in Alignments
http://rosalind.info/problems/osym/
"""
from __future__ import print_function
import argparse
import sys
import numpy

def fasta_iterator(infile):
    """
    Iterate over sequences in a FASTA file
    """
    started = False
    label = ""
    sequence = ""

    for line in infile:
        line = line.strip()
        if line[0] == '>':
            if started:
                yield label, sequence
            started = True
            label = line[1:]
            sequence = ""
        else:
            sequence += line

    if started:
        yield label, sequence

def global_match_alignment_score(seq1, seq2):
    """
    Find score of global alignment between seq1 + seq2
    """
    n = len(seq1)
    m = len(seq2)

    score = numpy.zeros((n+1, m+1, 4), dtype=numpy.int)
    match = 0
    s_ins = 1
    s_del = 2
    best = 3

    score[:, :, :] = 0
    for k in range(4):
        score[:, 0, k] = -1*numpy.arange(n+1)
        score[0, :, k] = -1*numpy.arange(m+1)

    for i in range(1, n+1):
        for j in range(1, m+1):
            matchscore = +1 if seq1[i-1] == seq2[j-1] else -1
            score[i, j, match] = score[i-1, j-1, best] + matchscore
            score[i, j, s_ins] = score[i-1, j, best] - 1
            score[i, j, s_del] = score[i, j-1, best] - 1
            score[i, j, best] = max(score[i, j, 0:3])

    return score

def isolation_matrix(seq1, seq2):
    """
    Find matrix of scores with individual base pairs
    forced into alignment
    """
    n = len(seq1)
    m = len(seq2)

    matchscore = numpy.zeros((n, m), dtype=numpy.int)
    for i in range(n):
        for j in range(m):
            matchscore[i, j] = +1 if seq1[i] == seq2[j] else -1

    # This is the best score so far involving a match at position i, j
    # which is the global alignment score for the left portion of
    # the strings up to and including positions i, j..
    left = global_match_alignment_score(seq1, seq2)[1:, 1:, 0]

    # This is the best score so far involving a match at position i, j
    # for the reverse strings
    # which is the global alignment score for the right portion of
    # the strings up to and including positions i, j..
    right = global_match_alignment_score(seq1[::-1], seq2[::-1])[1:, 1:, 0]

    # The score is the sum of the 
    score = left - matchscore
    for i in range(n):
        for j in range(m):
            score[i, j] = score[i, j] + right[n-1-i, m-1-j]

    return score

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for lab, seq in fasta_iterator(args.infile):
            sequences.append(seq)

        scores = isolation_matrix(sequences[1], sequences[0])
        print(numpy.max(scores))
        print(numpy.sum(scores))
