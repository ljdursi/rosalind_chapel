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


def failure_array(seq):
    n = len(seq)
    fa = numpy.zeros(n, dtype=numpy.int)

    prefixlen = 0
    prefixchar = seq[prefixlen]
    for i, base in zip(range(1, n), seq[1:]):
        while True:
            if base == prefixchar:
                prefixlen += 1
                prefixchar = seq[prefixlen]
                break
            if prefixlen == 0:
                break
            prefixlen = fa[prefixlen-1]
            prefixchar = seq[prefixlen]

        fa[i] = prefixlen

    return fa


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(f):
            fa = failure_array(sequence)
        
            print(" ".join([str(i) for i in fa]))
