#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import itertools

def count_chars(line):
    counts = defaultdict(int)
    for char in line:
        counts[char] += 1
    return counts

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        alphabet = lines[0].strip().split()
        k = int(lines[1].strip())

        for seq in itertools.product(alphabet, repeat=k):
            print("".join(seq))
