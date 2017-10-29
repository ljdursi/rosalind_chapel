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


def pdist(seq1, seq2):
    n = len(seq1)
    if len(seq2) != n:
        return -1
    
    ndiff = 0
    for b1, b2 in zip(seq1, seq2):
        if b1 != b2:
            ndiff += 1

    return 1.0 * ndiff / n


def distmatrix(sequences):
    nseq = len(sequences)
    d = numpy.zeros((nseq, nseq))

    for i in range(nseq):
        for j in range(i):
            d[i, j] = pdist(sequences[i], sequences[j])
            d[j, i] = d[i, j]

    return d


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(f):
            sequences.append(sequence)

        dists = distmatrix(sequences)
        for row in dists:
            print(" ".join([str(d) for d in row]))
