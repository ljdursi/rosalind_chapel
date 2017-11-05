#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import numpy


def scs(seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    # calculate length
    shortest = numpy.zeros((n+1, m+1), dtype=numpy.int)
    shortest[0, :] = range(m+1)
    shortest[:, 0] = range(n+1)

    for i in range(1, n+1):
        for j in range(1, m+1):
            if seq1[i-1] == seq2[j-1]:
                shortest[i, j] = shortest[i-1, j-1] + 1
            else:
                shortest[i, j] = min(shortest[i-1, j]+1, shortest[i, j-1]+1)

    # backtrack
    l = shortest[n, m]
    supersequence = [' ']*l
    i, j = n, m
    position = l

    while i > 0 or j > 0:
        up = shortest[i-1, j]
        left = shortest[i, j-1]
        diag = shortest[i-1, j-1]
        if diag == position-1 and i > 0 and j > 0 and seq1[i-1] == seq2[j-1]:
            position -= 1
            supersequence[position] = seq1[i-1]
            i -= 1
            j -= 1
        else:
            if left <= up and j > 0:
                position -= 1
                supersequence[position] = seq2[j-1]
                j -= 1
            else:
                position -= 1
                supersequence[position] = seq1[i-1]
                i -= 1

    return "".join(supersequence)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines]

        superseq = scs(lines[0], lines[1])
        print(superseq)
