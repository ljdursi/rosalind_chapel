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
    weights, aminos = zip(*rows)
    return aminos, weights


def nearest(numlist, query):
    n = len(numlist)

    loc = bisect.bisect_left(numlist, query)
    possibles = [loc-1, loc, loc+1]

    diffs = [(abs(query-numlist[i]), i) for i in possibles if i >= 0 and i < n]
    best = min(diffs)
    return best[1], best[0]


def infer_protein(prefixes, aminos, weights, tol=1.e-7):
    protein = ""
    startweight, endweight = prefixes[0], prefixes[-1]
    prefixes = prefixes[1:-1]
    n = len(prefixes)//2

    for i in range(n):
        best = (1000.0, None, None, None, None)

        for j, w in enumerate(weights):
            newstartweight = startweight + w
            newendweight = endweight - w

            test_prefix, deltap = nearest(prefixes, newstartweight)
            test_suffix, deltas = nearest(prefixes, newendweight)

            delta = max(deltap, deltas)
            if delta < best[0]:
                best = (delta, j, test_prefix, test_suffix, w)

            if delta < tol:
                break

        _, j, test_prefix, test_suffix, w = best
        protein = protein + aminos[j]
        startweight += w
        endweight -= w
        del prefixes[max(test_prefix, test_suffix)]
        del prefixes[min(test_prefix, test_suffix)]

    return protein


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-w', '--weightfile', default="weights.txt")
    parser.add_argument('-t', '--tolerance', type=float, default=1.e-7)
    args = parser.parse_args()

    aminos, weights = readtable(args.weightfile)
    with args.infile as f:
        lines = f.readlines()
        prefixes = [float(l.strip()) for l in lines]

        parentweight = prefixes.pop(0)
        prefixes.sort()

        result = infer_protein(prefixes, aminos, weights)
        print(result)
