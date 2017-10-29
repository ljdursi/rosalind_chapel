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


def composition(sequence, k):
    n = len(sequence)

    intseq = numpy.zeros(n, dtype=numpy.int)
    basevalue = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i, base in enumerate(sequence):
        intseq[i] = basevalue[base]

    kmer_value = 0
    for i in range(k):
        kmer_value *= 4
        kmer_value += intseq[i]

    nkmers = 4**k
    lastplace = 4**(k-1)
    counts = numpy.zeros(nkmers, dtype=numpy.int)
    for i in range(k, n):
        kmer_value -= lastplace*intseq[i-k]
        kmer_value = kmer_value*4 + intseq[i]
        counts[kmer_value] += 1

    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-k', '--kmersize', default=4, type=int)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(f):
            result = composition(sequence, args.kmersize)
            print(" ".join([str(i) for i in result]))
