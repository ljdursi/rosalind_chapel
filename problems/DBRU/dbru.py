#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys


def reverse_complement(line):
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rc = ''.join([complements[c] for c in reversed(line)])
    return rc


def readset(reads):
    rset = set([])
    for read in reads:
        rset.add(read)
        rset.add(reverse_complement(read))
    return rset


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]

        rset = readset(lines)
        for read in rset:
            print("({0}, {1})".format(read[:-1], read[1:]))
