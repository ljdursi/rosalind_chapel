#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import itertools

def factorial(n):
    if n <= 0:
        return 0

    result = 1
    for i in range(2,n+1):
        result *= i
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        line = "".join([l.strip() for l in lines])

        items = line.split()
        n = int(items[0])

        print(factorial(n))
        l = range(1,n+1)
        for perm in itertools.permutations(l):
            print(" ".join([str(p) for p in perm]))
