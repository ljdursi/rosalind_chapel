#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
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


def npk(n, k):
    result = 1
    for i in range(n-k+1, n+1):
        result *= i
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for label, sequence in fasta_iterator(f):
            counts = basecounts(sequence)

            aus = [counts['A'], counts['U']]
            cgs = [counts['C'], counts['G']]
            nmatches = npk(max(aus), min(aus))
            nmatches *= npk(max(cgs), min(cgs))
            print(nmatches)
