#!/usr/bin/env python
from __future__ import print_function
import argparse
import bisect
import sys
import numpy
from collections import defaultdict


def readtable(filename):
    aminos = []
    weights = []

    with open(filename, "r") as f:
        for line in f:
            items = line.strip().split()
            aminos.append(items[0])
            weights.append(float(items[1]))

    return aminos, weights


def nearest(numlist, query):
    n = len(numlist)

    loc = bisect.bisect_left(numlist, query)
    possibles = [loc-1, loc, loc+1]

    diffs = [(abs(query-numlist[i]), i) for i in possibles if i >= 0 and i < n]
    best = min(diffs)
    return best[1], best[0]


def spectrum_graph(peptides, aminos, weights, tol=1.e-6):
    naminos = len(aminos)
    peptides = sorted(peptides)
    graph = defaultdict(list)

    for i, peptide in enumerate(peptides):
        graph[i] = []
        for j in range(naminos):
            nweight = peptide + weights[j]
            neighbour, diff = nearest(peptides, nweight)
            if diff < tol:
                graph[i].append((aminos[j], neighbour))

    return graph


def longest_path(peptides, graph):
    n = len(peptides)

    best_pathlen = numpy.zeros(n, dtype=numpy.int)
    last = numpy.zeros(n, dtype=numpy.int) + -1
    longest = 0
    longestend = -1
    aminos = ['-']*n

    for i, peptide in enumerate(peptides):
        curbest = best_pathlen[i]
        for amino, neighbor in graph[i]:
            current = curbest + 1
            if current > best_pathlen[neighbor]:
                best_pathlen[neighbor] = current
                last[neighbor] = i
                aminos[neighbor] = amino
                if current > longest:
                    longest = current
                    longestend = neighbor

    longest_protein = [""]*longest
    idx, i = longest, longestend
    while i > -1 and idx > 0:
        idx -= 1
        longest_protein[idx] = aminos[i]
        i = last[i]

    result = "".join(longest_protein)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-w', '--weightfile', default="weights.txt")
    parser.add_argument('-t', '--tolerance', type=float, default=5.e-5)
    args = parser.parse_args()

    aminos, weights = readtable(args.weightfile)
    with args.infile as f:
        lines = f.readlines()
        peptides = [float(l) for l in lines]
        peptides.sort()

        graph = spectrum_graph(peptides, aminos, weights, args.tolerance)
        longest = longest_path(peptides, graph)
        print(longest)
