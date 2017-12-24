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


def readtable(filename):
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


def alignment(seq1, seq2, subscores, gap_parameter):
    n = len(seq1)
    m = len(seq2)

    score = numpy.zeros((n+1, m+1), dtype=numpy.int)
    move = numpy.zeros((n+1, m+1), dtype=numpy.int)
    move_diag = 3
    move_del = 2
    move_ins = 1
    move_start = 0

    move[:, :] = move_start

    for i in range(1, n+1):
        basei = seq1[i-1]
        for j in range(1, m+1):
            basej = seq2[j-1]

            match = score[i-1, j-1] + subscores[(basei, basej)]
            insertion = score[i-1, j] + gap_parameter
            deletion = score[i, j-1] + gap_parameter
            start = 0

            score[i, j], move[i, j] = max([(match, move_diag),
                                           (insertion, move_ins),
                                           (deletion, move_del),
                                           (start, move_start)])

    ij = numpy.argmax(score)
    i, j = numpy.unravel_index([ij], (n+1, m+1))
    i = i[0]
    j = j[0]
    maxscore = score[i, j]
    matched1 = ""
    matched2 = ""

    while True:
        if i <= 0 or j <= 0:
            break
        if move[i, j] == move_start:
            break
        if score[i, j] == 0:
            break

        if move[i, j] == move_diag:
            matched1 = seq1[i-1] + matched1
            matched2 = seq2[j-1] + matched2
            i -= 1
            j -= 1
        else:
            if move[i, j] == move_del:
                matched2 = seq2[j-1] + matched2
                j -= 1
            else:
                matched1 = seq1[i-1] + matched1
                i -= 1

    return maxscore, matched1, matched2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-s', '--scorefile', type=str, default="pam250.txt")
    parser.add_argument('-g', '--gapscore', type=int, default=-5)
    args = parser.parse_args()

    scoretable = readtable(args.scorefile)

    with args.infile as f:
        seqs = []
        for label, sequence in fasta_iterator(args.infile):
            seqs.append(sequence)

        score = alignment(seqs[0], seqs[1], scoretable, args.gapscore)
        print(score[0])
        print(score[1])
        print(score[2])
