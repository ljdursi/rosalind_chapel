#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import bisect


def readtable(filename="weights.txt"):
    rows = []

    with open(filename, "r") as f:
        for line in f:
            items = line.strip().split()
            amino = items[0]
            weight = float(items[1])
            rows.append((weight, amino))

    rows.sort()
    rows = [(0.0, None)] + rows + [(9999.0, None)]
    weights, aminos = zip(*rows)
    return aminos, weights


def infer_protein(prefixes, aminos, weights):
    prefixes.sort()

    protein = []
    for prefix1, prefix2 in zip(prefixes[:-1], prefixes[1:]):
        delta = prefix2-prefix1
        best = bisect.bisect_left(weights, delta)
        if abs(weights[best-1] - delta) < abs(weights[best] - delta):
            best = best - 1
        protein.append(aminos[best])

    return protein


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-w', '--weightfile', default="weights.txt")
    args = parser.parse_args()

    aminos, weights = readtable(args.weightfile)
    with args.infile as f:
        lines = f.readlines()
        prefixes = [float(l.strip()) for l in lines]

        result = infer_protein(prefixes, aminos, weights)
        print("".join(result))
