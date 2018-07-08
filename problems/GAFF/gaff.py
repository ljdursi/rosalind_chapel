#!/usr/bin/env python
"""
Global Alignment with Scoring Matrix and Affine Gap Penalty
http://rosalind.info/problems/gaff/

Useful background info: https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf
"""
from __future__ import print_function
import argparse
import sys
import numpy


def fasta_iterator(infile):
    """
    Iterate over sequences in a fasta
    """
    started = False
    label = ""
    sequence = ""

    for line in infile:
        line = line.strip()
        if line[0] == '>':
            if started:
                yield label, sequence
            started = True
            label = line[1:]
            sequence = ""
        else:
            sequence += line

    if started:
        yield label, sequence


def readtable(filename):
    """
    Read substitution score table,
    return as a dictionary over base pair tuples
    """
    table = {}
    with open(filename, "r") as f:
        header = f.readline()
        bases = header.split()
        for line in f:
            items = line.split()

            basei = items[0]
            vals = [int(i) for i in items[1:]]

            for basej, val in zip(bases, vals):
                table[(basei, basej)] = val

    return table


def alignment(seq1, seq2, subscores, gap, extend):
    """
    Global, Affine alignment
    """
    n = len(seq1)
    m = len(seq2)

    # Three score matrices:
    #     - M: best score for aligning seq1[1..i-1], seq2[1..j-1]
    #          with sequence ending in a match
    #     - X: best score with sequence ending in a gap in seq1
    #     - Y: best score with sequence ending in a gap in seq2
    # Three dimensions (e.g., M_ijk).  Third dimension:
    #    - i,j,0: score
    #    - i,j,1: i from previous move
    #    - i,j,2: j from previous move
    #    - i,j,3: matrix of previous move (0=m, 1=x, 2=y)
    scores = numpy.zeros((3, n+1, m+1, 4), dtype=numpy.int)
    m_mtx, x_mtx, y_mtx = 0, 1, 2
    m_score = scores[m_mtx, :, :, :]
    x_score = scores[x_mtx, :, :, :]
    y_score = scores[y_mtx, :, :, :]

    # boundary conditions:
    #   - no alignment of either length 0 that ends in a match,
    #     or gap in that sequence
    m_score[1:n+1, 0, 0] = numpy.iinfo(numpy.int).min/2
    m_score[0, 1:m+1, 0] = numpy.iinfo(numpy.int).min/2
    x_score[1:n+1, 0, 0] = gap + numpy.arange(1, n+1)*extend
    x_score[0, 1:m+1, 0] = numpy.iinfo(numpy.int).min/2
    y_score[1:n+1, 0, 0] = numpy.iinfo(numpy.int).min/2
    y_score[0, 1:m+1, 0] = gap + numpy.arange(1, m+1)*extend

    for i in range(1, n+1):
        basei = seq1[i-1]
        for j in range(1, m+1):
            basej = seq2[j-1]

            # match
            options = [(m_score[i-1, j-1, 0], i-1, j-1, m_mtx),
                       (x_score[i-1, j-1, 0], i-1, j-1, x_mtx),
                       (y_score[i-1, j-1, 0], i-1, j-1, y_mtx)]
            m_score[i, j, :] = max(options)
            m_score[i, j, 0] = m_score[i, j, 0] + subscores[(basei, basej)]

            # gap in seq 1
            options = [(m_score[i, j-1, 0] + gap + extend, i, j-1, m_mtx),
                       (x_score[i, j-1, 0] + extend, i, j-1, x_mtx),
                       (y_score[i, j-1, 0] + gap + extend, i, j-1, y_mtx)]
            x_score[i, j, :] = max(options)

            # gap in seq 1
            options = [(m_score[i-1, j, 0] + gap + extend, i-1, j, m_mtx),
                       (x_score[i-1, j, 0] + gap + extend, i-1, j, x_mtx),
                       (y_score[i-1, j, 0] + extend, i-1, j, y_mtx)]
            y_score[i, j, :] = max(options)

    # backtrack
    mtx = numpy.argmax(scores[:, n, m, 0])
    i = n
    j = m
    dist = scores[mtx, i, j, 0]
    augmented1 = ""
    augmented2 = ""

    while i > 0 or j > 0:
        move_i, move_j, move_mtx = scores[mtx, i, j, 1:]
        if mtx == m_mtx:
            augmented1 = seq1[i-1] + augmented1
            augmented2 = seq2[j-1] + augmented2
        elif mtx == x_mtx:
            augmented1 = '-' + augmented1
            augmented2 = seq2[j-1] + augmented2
        else:
            augmented1 = seq1[i-1] + augmented1
            augmented2 = '-' + augmented2

        mtx, i, j = move_mtx, move_i, move_j

    return dist, augmented1, augmented2

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-s', '--scorefile', type=str, default="blosum.txt")
    parser.add_argument('-g', '--gap_open', type=int, default=-10)
    parser.add_argument('-x', '--gap_extend', type=int, default=-1)
    args = parser.parse_args()

    scoretable = readtable(args.scorefile)

    with args.infile as f:
        seqs = []
        for lab, seq in fasta_iterator(args.infile):
            seqs.append(seq)

        score, aug1, aug2 = alignment(seqs[0], seqs[1], scoretable,
                                      args.gap_open, args.gap_extend)
        print(score)
        print(aug1)
        print(aug2)
