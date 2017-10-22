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
        lines = f.readlines()
        lines = [l.strip() for l in lines]

        nnodes = int(lines[0])
        nedges = len(lines[1:])

        print(nnodes-1-nedges)
