#!/usr/bin/env python
"""
Finding All Similar Motifs
http://rosalind.info/problems/ksim/
"""
from __future__ import print_function
import argparse
import sys
import numpy

def banded_alignment_mtx(score, seq1, seq2, bandwidth):
    """
    Return the score matrix of the banded edit alignment
    of seq1 and seq2
    """
    n = len(seq1)
    m = len(seq2)

    for i in range(1, n+1):
        start = max(1, i-bandwidth)
        finish = min(m, i+bandwidth)
        for j in range(start, finish+1):
            match = score[i-1, j-1] + 1 if seq1[i-1] != seq2[j-1] else score[i-1, j-1]
            insertion = score[i-1, j] + 1
            deletion = score[i, j-1] + 1

            score[i, j] = min(match, insertion, deletion)


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

        # set up scores matrix + boundary conditions once
        score = numpy.zeros((m+1, m+k+1), dtype=numpy.int)
        score[:, :] = 10*k
        score[0:k+1, 0] = range(0, k+1)
        score[0, 1:k+1] = range(1, k+1)

        for i in range(n-m+k+1):
            last = min(n, i+m+k)
            length = last-i+1
            banded_alignment_mtx(score, motif, sequence[i:last], k)

            for col in range(max(1, length-2*k-1), length):
                if score[m, col] <= k:
                    print(i+1, col)
