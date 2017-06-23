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

def signs_from_val(val, n):
    signs = []
    for _ in range(n):
        signs.append(2*(val % 2)-1)
        val /= 2
    return signs

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

        print(factorial(n) * 2**n)
        l = range(1,n+1)
        for perm in itertools.permutations(l):
            for signval in range(2**n):
                signs = signs_from_val(signval, n)
                print(" ".join([str(p*sign) for (p, sign) in zip(perm, signs)]))
