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
            q = numpy.array([float(item) for item in line.split()])
            p = 1. - q
            carrier = 2.*p*q
            print(" ".join([str(item) for item in carrier]))
