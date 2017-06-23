#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict

def hamming_distance(seq1, seq2):
    dist = 0
    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:
            dist += 1
    return dist

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        seq1 = lines[0].strip()
        seq2 = lines[1].strip()

        result = hamming_distance(seq1, seq2)
        print(result)
