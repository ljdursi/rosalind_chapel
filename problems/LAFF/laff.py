#!/usr/bin/env python
"""
Local Alignment with Scoring Matrix and Affine Gap Penalty
http://rosalind.info/problems/laff/

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
    #    - i,j,1: matrix of previous move (0=m, 1=x, 2=y)
    scores = numpy.zeros((n+1, m+1, 3, 2), dtype=numpy.int)
    m_mtx, x_mtx, y_mtx = 0, 1, 2
    m_score = scores[:, :, m_mtx, :]
    x_score = scores[:, :, x_mtx, :]
    y_score = scores[:, :, y_mtx, :]

    # boundary conditions:
    #   - no alignment of either length 0 that ends in a match,
    #     or gap in that sequence
    m_score[1:n+1, 0, 0] = numpy.iinfo(numpy.int32).min/2
    m_score[0, 1:m+1, 0] = numpy.iinfo(numpy.int32).min/2
    x_score[1:n+1, 0, 0] = 0
    x_score[0, 1:m+1, 0] = numpy.iinfo(numpy.int32).min/2
    y_score[1:n+1, 0, 0] = numpy.iinfo(numpy.int32).min/2
    y_score[0, 1:m+1, 0] = 0

    for i in range(1, n+1):
        basei = seq1[i-1]

        ## vectorize match and seq2 operations
        # match
        m_score[i, 1:m+1, 0] = m_score[i-1, 0:m, 0]
        m_score[i, 1:m+1, 1] = m_mtx

        x_better = numpy.where(x_score[i-1, 0:m, 0] > m_score[i, 1:m+1, 0])[0]
        m_score[i, x_better+1, 0] = x_score[i-1, x_better, 0]
        m_score[i, x_better+1, 1] = x_mtx

        y_better = numpy.where(y_score[i-1, 0:m, 0] > m_score[i, 1:m+1, 0])[0]
        m_score[i, y_better+1, 0] = y_score[i-1, y_better, 0]
        m_score[i, y_better+1, 1] = y_mtx

        # gap in seq 2
        y_score[i, 1:m+1, 0] = m_score[i-1, 1:m+1, 0] + gap + extend
        y_score[i, 1:m+1, 1] = m_mtx

        x_better = numpy.where((x_score[i-1, 1:m+1, 0] + gap + extend) > y_score[i, 1:m+1, 0])[0]
        y_score[i, x_better+1, 0] = x_score[i-1, x_better+1, 0] + gap + extend
        y_score[i, x_better+1, 1] = x_mtx

        y_better = numpy.where((y_score[i-1, 1:m+1, 0] + extend) > y_score[i, 1:m+1, 0])[0]
        y_score[i, y_better+1, 0] = y_score[i-1, y_better+1, 0] + extend
        y_score[i, y_better+1, 1] = y_mtx

        z_better = numpy.where(y_score[i, 1:m+1, 0] < 0)[0]
        y_score[i, z_better+1, 0] = 0
        y_score[i, z_better+1, 1] = 3

        for j in range(1, m+1):
            basej = seq2[j-1]

            # match
            m_score[i, j, 0] = m_score[i, j, 0] + subscores[(basei, basej)]
            if m_score[i, j, 0] < 0:
                m_score[i, j, 0] = 0
                m_score[i, j, 1] = 3

            # gap in seq 1
            x_score[i, j, 0] = m_score[i, j-1, 0]  + gap + extend
            x_score[i, j, 1] = m_mtx
            if x_score[i, j-1, 0] + extend > x_score[i, j, 0]:
                x_score[i, j, 0] = x_score[i, j-1, 0] + extend
                x_score[i, j, 1] = x_mtx
            if y_score[i, j-1, 0] + gap + extend > x_score[i, j, 0]:
                x_score[i, j, 0] = y_score[i, j-1, 0] + gap + extend
                x_score[i, j, 1] = y_mtx
            if 0 > x_score[i, j, 0]:
                x_score[i, j, 0] = 0
                x_score[i, j, 1] = 3

    # backtrack
    options = []
    for mtx in range(3):
        i, j = numpy.unravel_index(numpy.argmax(scores[:, :, mtx, 0]),
                                   scores[:, :, mtx, 0].shape)
        options.append((scores[i, j, mtx, 0], mtx, i, j))

    dist, mtx, i, j = max(options)
    sub1 = ""
    sub2 = ""
    sc = dist

    while i > 0 and j > 0 and sc > 0 and scores[i, j, mtx, 1] != 3:
        move_score, move_mtx = scores[i, j, mtx, :]
        if mtx == m_mtx:
            sub1 = seq1[i-1] + sub1
            sub2 = seq2[j-1] + sub2
            move_i, move_j = i-1, j-1
        elif mtx == x_mtx:
            sub2 = seq2[j-1] + sub2
            move_i, move_j = i, j-1
        else:
            sub1 = seq1[i-1] + sub1
            move_i, move_j = i-1, j

        sc, mtx, i, j = move_score, move_mtx, move_i, move_j

    return dist, sub1, sub2

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

        score, sub2, sub1 = alignment(seqs[1], seqs[0], scoretable,
                                      args.gap_open, args.gap_extend)
        print(score)
        print(sub1)
        print(sub2)
