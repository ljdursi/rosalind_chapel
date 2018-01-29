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
    with open(filename, "r") as table_f:
        header = table_f.readline()
        bases = header.split()
        for line in table_f:
            items = line.split()

            basei = items[0]
            vals = [int(i) for i in items[1:]]

            for basej, val in zip(bases, vals):
                table[(basei, basej)] = val

    return table


def alignment(seq1, seq2, subs, gap_parameter):
    n = len(seq1)
    m = len(seq2)

    match_scores = numpy.zeros((n+1, m+1), dtype=numpy.int)
    ins_scores = numpy.zeros((n+1, m+1), dtype=numpy.int)
    del_scores = numpy.zeros((n+1, m+1), dtype=numpy.int)

    minint = numpy.iinfo(numpy.int32).min

    match_scores[0, 1:m+1] = minint
    match_scores[1:n+1, 0] = minint
    del_scores[1:n+1, 0] = minint
    del_scores[0, 1:m+1] = gap_parameter
    ins_scores[1:n+1, 0] = gap_parameter
    ins_scores[0, 1:m+1] = minint

    for i in range(1, n+1):
        bi = seq1[i-1]
        for j in range(1, m+1):
            bj = seq2[j-1]

            match_scores[i, j] = max([match_scores[i-1, j-1],
                                      ins_scores[i-1, j-1],
                                      del_scores[i-1, j-1]]) + subs[(bi, bj)]

            ins_scores[i, j] = max([match_scores[i-1, j] + gap_parameter,
                                    del_scores[i-1, j] + gap_parameter,
                                    ins_scores[i-1, j]])

            del_scores[i, j] = max([match_scores[i, j-1] + gap_parameter,
                                    ins_scores[i, j-1] + gap_parameter,
                                    del_scores[i, j-1]])

    return max([match_scores[n, m], ins_scores[n, m], del_scores[n, m]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-s', '--scorefile', type=str, default="blosum.txt")
    parser.add_argument('-g', '--gapscore', type=int, default=-5)
    args = parser.parse_args()

    scoretable = readtable(args.scorefile)

    with args.infile as f:
        seqs = []
        for seq_label, seq in fasta_iterator(args.infile):
            seqs.append(seq)

        result = alignment(seqs[0], seqs[1], scoretable, args.gapscore)
        print(result)
