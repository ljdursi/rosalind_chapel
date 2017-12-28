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


def motif_alignment(motif, seq2):
    n = len(motif)
    m = len(seq2)

    score = numpy.zeros((n+1, m+1), dtype=numpy.int)
    move = numpy.zeros((n+1, m+1), dtype=numpy.int)
    move_diag = 3
    move_del = 2
    move_ins = 1

    score[:, :] = 0
    score[:, 0] = -1*numpy.arange(n+1)
    move[:, :] = 0

    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score[i-1, j-1] + 1
            mismtch = score[i-1, j-1] - 1
            insertion = score[i-1, j] - 1
            deletion = score[i, j-1] - 1
            options = [(insertion, move_ins), (deletion, move_del)]
            if motif[i-1] == seq2[j-1]:
                score[i, j], move[i, j] = max([(match, move_diag)] + options)
            else:
                score[i, j], move[i, j] = max([(mismtch, move_diag)] + options)

    # backtrack
    j = numpy.argmax(score[n, :])
    dist = score[n, j]
    augmented1 = ""
    augmented2 = ""
    i = n

    while i > 0 and j > 0:
        if move[i, j] == move_diag:
            augmented1 = motif[i-1] + augmented1
            augmented2 = seq2[j-1] + augmented2
            i -= 1
            j -= 1
        else:
            if move[i, j] == move_del:
                augmented2 = seq2[j-1] + augmented2
                augmented1 = '-' + augmented1
                j -= 1
            else:
                augmented1 = motif[i-1] + augmented1
                augmented2 = '-' + augmented2
                i -= 1

    return dist, augmented1, augmented2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(args.infile):
            sequences.append(sequence)

        d, aug1, aug2 = motif_alignment(sequences[1], sequences[0])
        print(d)
        print(aug2)
        print(aug1)
