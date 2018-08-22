#!/usr/bin/env python
"""
Maximum gap alignment
http://rosalind.info/problems/mgap/
"""
from __future__ import print_function
import argparse
import sys
import numpy


def fasta_iterator(infile):
    """
    Simple parser of fasta files
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


def edit_dist(seq1, seq2, match_score=+10, mismatch_score=-3, gap_score=-1):
    n = len(seq1)
    m = len(seq2)

    score = numpy.zeros((n+1, m+1), dtype=numpy.int)
    ngaps = numpy.zeros((n+1, m+1), dtype=numpy.int)
    score[1:n+1, 0] = numpy.arange(1, n+1)*gap_score
    ngaps[1:n+1, 0] = numpy.arange(1, n+1)
    score[0, 1:m+1] = numpy.arange(1, m+1)*gap_score
    ngaps[0, 1:m+1] = numpy.arange(1, m+1)
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score[i-1, j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            insertion = score[i-1, j] + gap_score
            deletion = score[i, j-1] + gap_score
            (score[i, j], ngaps[i, j]) = max((match, ngaps[i-1, j-1]),
                                             (insertion, ngaps[i-1, j]+1),
                                             (deletion, ngaps[i, j-1]+1))

    return score[n, m], ngaps[n, m]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('--match', '-m', type=int, default=+10)
    parser.add_argument('--mismatch', '-d', type=int, default=-3)
    parser.add_argument('--gap', '-g', type=int, default=-1)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(args.infile):
            sequences.append(sequence)

        bestscore, bestgaps = edit_dist(sequences[0], sequences[1],
                                        match_score=args.match, mismatch_score=args.mismatch,
                                        gap_score=args.gap)
        print(bestgaps)
