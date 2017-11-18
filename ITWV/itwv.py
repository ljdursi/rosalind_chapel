#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import numpy


baseidx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def is_interleaving(superseq, seq1, seq2):
    n = len(seq1)
    m = len(seq2)
    l = len(superseq)
    if not n+m == l:
        return False

    scores = numpy.zeros((n+1, m+1), dtype=numpy.int)
    for i in range(1, n+1):
        if seq1[i-1] == superseq[i-1]:
            scores[i, 0] = scores[i-1, 0] + 1

    for j in range(1, m+1):
        if seq2[j-1] == superseq[j-1]:
            scores[0, j] = scores[0, j-1] + 1

    for i in range(1, n+1):
        for j in range(1, m+1):
            idx = i + j - 1
            left = scores[i, j-1]
            if seq2[j-1] == superseq[idx]:
                left += 1
            down = scores[i-1, j]
            if seq1[i-1] == superseq[idx]:
                down += 1
            scores[i, j] = max(left, down)

    return(scores[n, m] == l)


def calc_counts(sequence):
    n = len(sequence)
    nb = len(baseidx)
    counts = numpy.zeros((nb, n+1), dtype=numpy.int)
    for i, base in enumerate(sequence):
        counts[:, i+1] = counts[:, i]
        counts[baseidx[base]][i+1] = counts[baseidx[base]][i] + 1
    return counts


def find_windows(sseq_counts, sseq_len, seq1, seq2):
    totlen = len(seq1) + len(seq2)
    basecounts = numpy.zeros(len(baseidx), dtype=numpy.int)
    for base in seq1 + seq2:
        basecounts[baseidx[base]] += 1

    windows = []
    for idx in range(0, sseq_len-totlen+1):
        delta = sseq_counts[:, idx+totlen] - sseq_counts[:, idx]
        if numpy.all(basecounts == delta):
            windows.append(idx)
    return windows


def contains_interleaving(counts, sseq, seq1, seq2):
    windows = find_windows(counts, len(sseq), seq1, seq2)

    totlen = len(seq1) + len(seq2)
    found = False
    for i in windows:
        found = is_interleaving(sseq[i:i+totlen], seq1, seq2)
        if found:
            break

    return found


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = [line.strip() for line in f.readlines()]
        superseq = lines[0]
        counts = calc_counts(superseq)
        subseqs = lines[1:]

        n = len(subseqs)
        mtx = numpy.zeros((n, n), dtype=numpy.int)
        for i in range(n):
            for j in range(i+1):
                result = contains_interleaving(counts, superseq,
                                               subseqs[i], subseqs[j])
                if result:
                    mtx[i, j] = 1
                    mtx[j, i] = 1

        for i in range(n):
            print(" ".join([str(mj) for mj in mtx[i, :]]))
