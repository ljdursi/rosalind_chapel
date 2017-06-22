#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict

def generalized_fib(ngen, off_per_gen):
    last, current = 0, 1
    for i in range(ngen-1):
        last, current = current, current + off_per_gen*last
    return current

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        line = "".join([l.strip() for l in lines])

        items = line.split()
        n, k = int(items[0]), int(items[1])

        result = generalized_fib(n, k)
        print(result)
