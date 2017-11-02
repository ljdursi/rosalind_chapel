#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import numpy


def apply_permutation(permutation, series):
    n = series.size
    result = numpy.zeros(n, dtype=numpy.int)

    for i in range(n):
        result[i] = permutation[series[i]-1]
    return result


def inverse_permutation(permutation):
    n = permutation.size
    result = numpy.zeros(n, dtype=numpy.int)

    for i in range(n):
        result[permutation[i]-1] = i + 1
    return result


def adjacent(i, j):
    return i+1 == j or j+1 == i


def breakpoints(sequence):
    n = sequence.size
    bps = []
    for i in range(1, n):
        if not adjacent(sequence[i-1], sequence[i]):
            bps.append(i)
    return bps


def potential_reversals(sequence, bkpts):
    max_new_adj = -1
    nbps = len(bkpts)
    potentials = []

    for i in range(1, nbps):
        for j in range(0, i):
            if bkpts[i] == bkpts[j] + 1:
                continue
            new_l_adj = adjacent(sequence[bkpts[j]-1], sequence[bkpts[i]-1])
            new_r_adj = adjacent(sequence[bkpts[j]], sequence[bkpts[i]])
            n_new_adj = (1 if new_l_adj else 0) + (1 if new_r_adj else 0)
            if n_new_adj > max_new_adj:
                potentials = []
                max_new_adj = n_new_adj
            if n_new_adj == max_new_adj:
                potentials.append((j, i))

    return potentials


def exhaustive_reversal_distance(sequence, best=None):
    bkpts = breakpoints(sequence)
    n = len(bkpts)
    if n == 0:
        return 0

    if not best:
        best = 2*n

    if best == 0:
        return 2*n

    moves = potential_reversals(sequence, bkpts)
    for leftbkpt, rightbkpt in moves:
        newsequence = numpy.copy(sequence)
        start, end = bkpts[leftbkpt], bkpts[rightbkpt]
        newsequence[start:end] = newsequence[start:end][::-1]
        dist = exhaustive_reversal_distance(newsequence, best-1)
        if dist < best:
            best = dist
        if best <= (n+1)//2:
            break

    return 1 + best


def reversal_distance(perm1, perm2):
    n = perm2.size
    invp1 = inverse_permutation(perm1)
    normalized = apply_permutation(invp1, perm2)
    extended = numpy.zeros(normalized.size+2, dtype=numpy.int)
    extended[1:n+1] = normalized
    extended[n+1] = n+1

    return exhaustive_reversal_distance(extended)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines]
        lines = [l for l in lines if len(l) != 0]

        rev_distances = []
        for line1, line2 in zip(lines[::2], lines[1::2]):
            perm1 = numpy.array([int(i) for i in line1.split()])
            perm2 = numpy.array([int(i) for i in line2.split()])
            rev_distances.append(reversal_distance(perm1, perm2))

        print(" ".join([str(i) for i in rev_distances]))
