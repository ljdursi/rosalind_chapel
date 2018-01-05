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
    wobblematchings = {"G": "U", "U": "G"}
    base1 = sequence[-1]
    base2 = matchings[base1]
    can_wobble = base1 in wobblematchings
    if can_wobble:
        base3 = wobblematchings[base1]

    return [i for i in range(len(sequence)-4)
            if sequence[i] == base2 or (can_wobble and sequence[i] == base3)]


def non_crossing_matchings(seq):
    m = n = len(seq)
    matchings = numpy.zeros((m+1, n+1), dtype='object')
    matchings[:, :] = 0
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
                nmatch += (nleft * nright)

            matchings[length, start] = nmatch

    return matchings[m, 0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            sequence = line.strip()
            print(non_crossing_matchings(sequence))
