#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import numpy
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-k', '--overlap', type=int, default=3)
    args = parser.parse_args()

    with args.infile as f:
        prefixes = defaultdict(list)
        reads = []
        for label, sequence in fasta_iterator(f):
            reads.append((label,sequence))
            prefixes[sequence[0:args.overlap]].append((label,sequence))

        for label, sequence in reads:
            suffix = sequence[-args.overlap:]
            for neighbor_label, neighbor_sequence in prefixes[suffix]:
                if neighbor_label != label:
                    print(label, neighbor_label)
