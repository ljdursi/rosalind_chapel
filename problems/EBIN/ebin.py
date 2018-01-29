#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = [line.strip() for line in f.readlines()]
        npop = int(lines[0])
        probs = [float(item) for item in lines[1].split()]
        print(" ".join([str(npop*p) for p in probs]))
