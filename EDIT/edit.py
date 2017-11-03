#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import numpy


def fasta_iterator(infile):
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


def edit_dist(seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    score = numpy.zeros((n+1, m+1), dtype=numpy.int)
    score[1:n+1, 0] = range(1, n+1)
    score[0, 1:m+1] = range(1, m+1)
    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score[i-1, j-1]
            mismatch = score[i-1, j-1] + 1
            insertion = score[i-1, j] + 1
            deletion = score[i, j-1] + 1
            if seq1[i-1] == seq2[j-1]:
                score[i, j] = min(match, insertion, deletion)
            else:
                score[i, j] = min(mismatch, insertion, deletion)

    return score[n, m]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(args.infile):
            sequences.append(sequence)

        result = edit_dist(sequences[0], sequences[1])
        print(result)
