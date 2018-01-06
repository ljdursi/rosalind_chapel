#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys


def double_factorial_mod(n, mod):
    result = 1
    for i in range(1, n, 2):
        result *= i
        result %= mod
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-m', '--modulo', type=int, default=1000000)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            n = int(line)
            print(double_factorial_mod(2*n-1, args.modulo))
