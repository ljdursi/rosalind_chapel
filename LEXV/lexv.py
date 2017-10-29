#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys


def substrings(alphabet, k, first=True):
    if k == 0:
        yield ""
    else:
        if not first:
            yield ""
        for head in alphabet:
            tails = substrings(alphabet, k-1, False)
            for tail in tails:
                yield head + tail


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        alphabet = lines[0].strip().split()
        k = int(lines[1].strip())

        for seq in substrings(alphabet, k):
            print(seq)
