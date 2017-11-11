#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys


def n_choose_k_mod_m(n, k, modulus):
    if k < 0 or k > n:
        return 0

    if n == k or k == 0:
        return 1

    num, den = 1, 1
    for i in range(1, min(k, n-k) + 1):
        num *= n
        den *= i
        n -= 1
    return (num // den) % modulus


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-m', '--modulus', type=int, default=1000000)
    parser.add_argument('-n', type=int, default=2)
    args = parser.parse_args()
    modulus = args.modulus

    with args.infile as f:
        lines = f.readlines()
        line = "".join([l.strip() for l in lines])

        items = line.split()
        n, m = int(items[0]), int(items[1])

        res = 0
        for k in range(m, n+1):
            res = (res + n_choose_k_mod_m(n, k, modulus)) % modulus

        print(res)
