#!/usr/bin/env python
"""
Finding All Similar Motifs
http://rosalind.info/problems/ksim/
"""
from __future__ import print_function
import argparse
import sys
import numpy


def banded_alignment_mtx(seq1, seq2, bandwidth):
    """
    Return the score matrix of the banded edit alignment
    of seq1 and seq2
    """
    n = len(seq1)
    m = len(seq2)

    score = numpy.zeros((n+1, m+1), dtype=numpy.int)
    score[:, :] = 10*bandwidth
    score[0:bandwidth+1, 0] = range(0, bandwidth+1)
    score[0, 1:bandwidth+1] = range(1, bandwidth+1)

    for i in range(1, n+1):
        start = max(1, i-bandwidth)
        finish = min(m, i+bandwidth)
        for j in range(start, finish+1):
            match = score[i-1, j-1] + 1 if seq1[i-1] != seq2[j-1] else score[i-1, j-1]
            insertion = score[i-1, j] + 1
            deletion = score[i, j-1] + 1

            score[i, j] = min(match, insertion, deletion)

    return score


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = []
        for line in f:
            lines.append(line.strip())

        k, motif, sequence = int(lines[0]), lines[1], lines[2]
        m = len(motif)
        n = len(sequence)

        for i in range(n-m+k+1):
            if i % 100 == 0:
                print(i, file=sys.stderr)

            last = min(n, i+m+k-1)
            # length = last-i+1
            scores = banded_alignment_mtx(motif, sequence[i:last+1], k)
            _, length = scores.shape

            for col in range(max(1, length-2*k-1), length):
                if scores[m, col] <= k:
                    print(i+1, col)
