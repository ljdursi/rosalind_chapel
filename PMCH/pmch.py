#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import math
from collections import defaultdict


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


def basecounts(seq):
    counts = defaultdict(int)
    for base in seq:
        counts[base] += 1

    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for label, sequence in fasta_iterator(f):
            counts = basecounts(sequence)
            assert counts['A'] == counts['U']
            assert counts['C'] == counts['G']

            nmatches = math.factorial(counts['A'])
            nmatches *= math.factorial(counts['C'])
            print(nmatches)
