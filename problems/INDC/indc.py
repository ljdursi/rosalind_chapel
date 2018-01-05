#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import math
from scipy.stats import binom


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-n', '--nalleles', type=int, default=2)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            line = line.strip()
            items = line.split()

            n = int(items[0])

            dist = binom(2*n, 0.5)
            bsf = [math.log10(binom.sf(i, 2*n, 0.5)) for i in range(2*n)]
            print(" ".join([str(r) for r in bsf]))
