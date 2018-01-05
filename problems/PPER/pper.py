#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys


def partial_factorial_modulus(n, k, modulus):
    if n <= 0:
        return 0

    result = 1
    for i in range(n, n-k, -1):
        result *= i
        result %= modulus
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-m', '--modulus', type=int, default=1000000)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        line = "".join([l.strip() for l in lines])

        items = line.split()
        n = int(items[0])
        k = int(items[1])

        print(partial_factorial_modulus(n, k, args.modulus))
