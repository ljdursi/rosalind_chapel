#!/usr/bin/env python
from __future__ import print_function
import argparse
import bisect
import sys
import numpy


def readtable(filename):
    amino_weights = {}

    with open(filename, "r") as f:
        for line in f:
            items = line.strip().split()
            amino = items[0]
            weight = float(items[1])
            amino_weights[amino] = weight

    return amino_weights


def nearest(numlist, query):
    n = len(numlist)

    loc = bisect.bisect_left(numlist, query)
    possibles = [loc-1, loc, loc+1]

    diffs = [(abs(query-numlist[i]), i) for i in possibles if i >= 0 and i < n]
    best = min(diffs)
    return best[1], best[0]


def differences(a, b):
    n = len(a)
    m = len(b)
    diffs = [0.0]*(n*m)

    i = 0
    for aa in a:
        for bb in b:
            diffs[i] = aa-bb
            i += 1

    diffs.sort()
    return diffs


def multiplicities(l, tol):
    last = -999999999999
    count = 0
    multiplicities = []
    for item in l:
        if abs(item-last) < tol:
            count += 1
        else:
            multiplicities.append((count, last))
            last = item
            count = 1

    return multiplicities[1:] + [(count, last)]


def spectrum(protein, amino_weights):
    weights = numpy.array([amino_weights[aa] for aa in protein])
    total_weight = numpy.sum(weights)

    prefixes = numpy.cumsum(weights)
    suffixes = total_weight - prefixes

    spec = [total_weight] + list(prefixes) + list(suffixes[:-1])
    spec.sort()
    return spec


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-w', '--weightfile', default="weights.txt")
    parser.add_argument('-t', '--tolerance', type=float, default=1.e-6)
    args = parser.parse_args()

    amino_weights = readtable(args.weightfile)
    with args.infile as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines]

        n = int(lines[0])
        proteins = lines[1:n+1]
        peaks = [float(l) for l in lines[n+1:]]
        peaks.sort()

        best = (-1, None)
        for protein in proteins:
            spec = spectrum(protein, amino_weights)
            diffs = differences(spec, peaks)
            multiplicity = multiplicities(diffs, args.tolerance)
            maxmult = max(multiplicity)
            if (maxmult[0], protein) > best:
                best = (maxmult[0], protein)

        print(best[0])
        print(best[1])
