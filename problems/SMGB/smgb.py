#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import numpy


def fasta_iterator(infile):
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


def semiglobal_alignment(seq1, seq2, match, mismatch, gap):
    n = len(seq1)
    m = len(seq2)

    score = numpy.zeros((n+1, m+1), dtype=numpy.int)
    score[0, :] = 0
    score[:, 0] = 0

    move = numpy.zeros((n+1, m+1), dtype=numpy.int)
    move_diag = 4
    move_del = 2
    move_ins = 1

    move[:, :] = move_diag
    move[0, :] = move_del
    move[:, 0] = move_ins

    for i in range(1, n+1):
        basei = seq1[i-1]
        del_gap = gap if i < n else 0

        for j in range(1, m+1):
            ins_gap = gap if j < m else 0

            match_score = score[i-1, j-1]
            if basei == seq2[j-1]:
                match_score += match
            else:
                match_score += mismatch

            insertion_score = score[i-1, j] + ins_gap
            deletion_score = score[i, j-1] + del_gap
            score[i, j], move[i, j] = max([(insertion_score, move_ins),
                                           (deletion_score, move_del),
                                           (match_score, move_diag)])

    i = n
    j = m
    maxscore = score[i, j]
    matched1 = ""
    matched2 = ""

    while True:
        if i <= 0 or j <= 0:
            break

        if move[i, j] == move_diag:
            matched1 = seq1[i-1] + matched1
            matched2 = seq2[j-1] + matched2
            i -= 1
            j -= 1
        else:
            if move[i, j] == move_del:
                matched2 = seq2[j-1] + matched2
                matched1 = "-" + matched1
                j -= 1
            else:
                matched1 = seq1[i-1] + matched1
                matched2 = "-" + matched2
                i -= 1

    return maxscore, matched1, matched2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-g', '--gapscore', type=int, default=-1)
    parser.add_argument('-m', '--matchscore', type=int, default=1)
    parser.add_argument('-M', '--mismatchscore', type=int, default=-1)
    args = parser.parse_args()

    with args.infile as f:
        seqs = []
        for label, sequence in fasta_iterator(args.infile):
            seqs.append(sequence)

        score = semiglobal_alignment(seqs[1], seqs[0], args.matchscore,
                                  args.mismatchscore, args.gapscore)
        print(score[0])
        print(score[2])
        print(score[1])
