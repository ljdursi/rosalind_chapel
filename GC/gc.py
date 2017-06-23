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

def sequence_gc(sequence):
    n = len(sequence)
    gc = 0
    for base in sequence:
        if base in ['G', 'C']:
            gc += 1
    return 100.*gc/n

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        maxgc = -1
        maxlabel = "none"
        for label, sequence in fasta_iterator(f):
            gc = sequence_gc(sequence)
            if gc > maxgc:
                maxgc, maxlabel = gc, label
        print(maxlabel)
        print(maxgc)
