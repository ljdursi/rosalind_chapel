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


def delta_base_counts(sequence, base1, base2):
    n = len(sequence)
    deltas = numpy.zeros(n, dtype=numpy.int)

    def inc(char, c1, c2):
        if char == c1:
            return +1
        elif char == c2:
            return -1
        else:
            return 0

    deltas[0] = inc(sequence[0], base1, base2)
    for i in range(1, n):
        base = sequence[i]
        deltas[i] = deltas[i-1] + inc(base, base1, base2)

    return deltas


def possible_partitions(sequence, a_minus_u, c_minus_g):
    matchings = {"A": "U", "C": "G", "G": "C", "U": "A"}
    base1 = sequence[0]
    base2 = matchings[base1]

    partitions = []
    for i in range(1, len(sequence)):
        if sequence[i] == base2:
            if a_minus_u[i] == 0 and c_minus_g[i] == 0:
                partitions.append(i)

    return partitions


def non_crossing_matchings(seq, modulus=1000000):
    a_minus_u = delta_base_counts(seq, "A", "U")
    c_minus_g = delta_base_counts(seq, "C", "G")

    m = n = len(seq)
    matchings = numpy.zeros((m+1, n+1), dtype=numpy.int)
    matchings[0, :] = 1

    for length in range(2, m+1, 2):
        for start in range(n-length+1):
            part_seq = seq[start:start+length]
            part_amu = a_minus_u[start:start+length].copy()
            part_cmg = c_minus_g[start:start+length].copy()
            if start > 0:
                part_amu -= a_minus_u[start-1]
                part_cmg -= c_minus_g[start-1]

            partitions = possible_partitions(part_seq, part_amu, part_cmg)
            nmatch = 0
            for p in partitions:
                left_start, left_len = start + 1, p - 1
                right_start, right_len = start + p + 1, length - p - 1
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
