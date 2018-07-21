#!/usr/bin/env python
"""
Multiple Alignment
http://rosalind.info/problems/mult/
"""
from __future__ import print_function
import argparse
import sys
import itertools
import numpy

def fasta_iterator(infile):
    """
    Iterate over sequences in a fasta file
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


def delta_score(align_bases):
    """
    Given a column of the alignment, eg ['A','A','T','-'],
    return the change in score - -1*number of pairwise mismatches
    """
    score = 0
    for base1, base2 in itertools.combinations(align_bases, 2):
        if base1 != base2:
            score -= 1
    return score


def delta_idxs_from_move(move, dimension):
    partial_moves = 2**numpy.arange(dimension)
    pmove = -move
    delta_idxs = [0]*dimension
    for i in range(dimension-1,-1,-1):
        if pmove // partial_moves[i] == 1:
            pmove -= partial_moves[i]
            delta_idxs[i] = -1
    return delta_idxs


def multiple_alignment(seqs):
    dimension = len(seqs)
    sizes = tuple(len(seq) for seq in seqs)
    mtx_sizes = tuple([s+1 for s in sizes])

    score = numpy.zeros(mtx_sizes, dtype=numpy.int)
    move = numpy.zeros(mtx_sizes, dtype=numpy.int)

    partial_moves = 2**numpy.arange(dimension)

    for idx, _ in numpy.ndenumerate(score):
        if idx == (0, 0, 0, 0):
            continue
        score[idx] = -999999
        for delta in itertools.product([-1, 0], repeat=dimension):
            prev_idx = tuple(i + d for i, d in zip(idx, delta))
            if any(i < 0 for i in prev_idx):
                continue
            align_bases = tuple('-' if delta[d] == 0 else seqs[d][idx[d]-1] for d in range(dimension))
            trial_score = score[prev_idx] + delta_score(align_bases)
            trial_move = sum([partial_moves[i]*delta[i] for i in range(dimension)])

            if trial_score > score[idx]:
                score[idx] = trial_score
                move[idx] = trial_move

    # backtrack
    idx = sizes
    dist = score[idx]
    augmenteds = [""]*dimension

    while not all(i == 0 for i in idx):
        delta = delta_idxs_from_move(move[idx], dimension)
        prev_idx = tuple(i + d for i, d in zip(idx, delta))
        align_bases = tuple('-' if delta[d] == 0 else seqs[d][idx[d]-1] for d in range(dimension))

        augmenteds = [align + aug for (align, aug) in zip(align_bases, augmenteds)]
        idx = prev_idx

    return dist, augmenteds


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(args.infile):
            sequences.append(sequence)

        d, augs = multiple_alignment(sequences)
        print(d)
        for aug in augs:
            print(aug)
