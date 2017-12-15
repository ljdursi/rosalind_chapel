#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import numpy

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            a = numpy.array([float(item) for item in line.split()])
            q = numpy.sqrt(a)
            b = 1.0-(1.0-q)**2
            print(" ".join([str(item) for item in b]))
