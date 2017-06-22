#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict

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
        line = "".join([l.strip() for l in lines])

        result = count_chars(line)
        print(result['A'], result['C'], result['G'], result['T'])
