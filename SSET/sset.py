#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys


def n_to_the_k_mod(n, k, modulus):
    if k <= 0:
        return 1

    result = 1
    for i in range(k):
        result *= n
        result %= modulus
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-m', '--modulus', type=int, default=1000000)
    parser.add_argument('-n', type=int, default=2)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        line = "".join([l.strip() for l in lines])

        items = line.split()
        k = int(items[0])

        print(n_to_the_k_mod(args.n, k, args.modulus))
