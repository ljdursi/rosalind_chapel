#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import numpy


def lcs(seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    # calculate length
    longest = numpy.zeros((n+1, m+1), dtype=numpy.int)
    for i in range(1, n+1):
        for j in range(1, m+1):
            if seq1[i-1] == seq2[j-1]:
                longest[i, j] = longest[i-1, j-1] + 1
            else:
                longest[i, j] = max(longest[i-1, j], longest[i, j-1])

    # backtrack
    l = longest[n, m]

    subsequence = [0]*l
    i, j = n, m
    position = l

    while i > 0 and j > 0:
        up = longest[i-1, j]
        left = longest[i, j-1]
        if up < position and left < position:
            position -= 1
            subsequence[position] = seq1[i-1]
            i -= 1
            j -= 1
        else:
            if left > up:
                j -= 1
            else:
                i -= 1

    return subsequence


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines]

        sequence = [int(item) for item in lines[1].split()]
        nitems = len(sequence)

        increasing = lcs(sequence, range(1, nitems+1))
        print(" ".join([str(i) for i in increasing]))

        decreasing = lcs(sequence, range(nitems, 0, -1))
        print(" ".join([str(i) for i in decreasing]))
