#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict

def reversecomplement(line):
    complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    rc = ''.join([complements[c] for c in reversed(line)])
    return rc

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        line = "".join([l.strip() for l in lines])

        result = reversecomplement(line)
        print(result)
