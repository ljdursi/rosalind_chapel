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
        line = "".join([l.strip() for l in lines])

        items = line.split()
        k, m, n = int(items[0]), int(items[1]), int(items[2])
        tot = k + m + n

        p = (1./tot)*(k + m/(tot-1.0)*(k + 0.75*(m-1.0) + 0.5*n) + n/(tot-1.0)*(k + 0.5*m))
        print(p)
