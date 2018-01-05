#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys


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


def subsequence_indices(sequence, subsequence):
    seq_cursor = 0
    subseq_positions = []

    for i, base in enumerate(subsequence): 
        while not sequence[seq_cursor] == base:
            seq_cursor += 1
        subseq_positions.append(seq_cursor)
        seq_cursor += 1

    return subseq_positions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(f):
            sequences.append(sequence)
        
        indices = subsequence_indices(sequences[0], sequences[1])
        print(" ".join([str(i+1) for i in indices]))
