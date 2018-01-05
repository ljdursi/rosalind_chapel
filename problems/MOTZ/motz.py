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


def possible_partitions(sequence):
    matchings = {"A": "U", "C": "G", "G": "C", "U": "A"}
    base1 = sequence[-1]
    base2 = matchings[base1]

    return [i for i in range(len(sequence)-1) if sequence[i] == base2]


def non_crossing_matchings(seq, modulus=1000000):
    m = n = len(seq)
    matchings = numpy.zeros((m+1, n+1), dtype=numpy.int)
    matchings[0, :] = 1
    matchings[1, :] = 1

    for length in range(2, m+1):
        for start in range(n-length+1):
            nmatch = matchings[length-1, start]
            partitions = possible_partitions(seq[start:start+length])
            for p in partitions:
                left_start, left_len = start, p
                right_start, right_len = start + p + 1, length - p - 2
                nleft = matchings[left_len, left_start]
                nright = matchings[right_len, right_start]
                nmatch += (nleft * nright) % modulus

            matchings[length, start] = nmatch % modulus

    return matchings[m, 0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-m', '--modulus', default=1000000, type=int)
    args = parser.parse_args()

    with args.infile as f:
        for label, sequence in fasta_iterator(args.infile):
            print(non_crossing_matchings(sequence, args.modulus))
