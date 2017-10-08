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


def transition_transversion(seq1, seq2):
    basetype = {"A": "purine", "G": "purine",
                "C": "pyramidine", "T": "pyramidine"}

    transitions = 0
    transversions = 0
    for base1, base2 in zip(seq1, seq2):
        if not base1 == base2:
            if basetype[base1] == basetype[base2]:
                transitions += 1
            else:
                transversions += 1

    return transitions, transversions


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        seqs = []
        for label, sequence in fasta_iterator(f):
            seqs.append(sequence)
        
        transition, transversion = transition_transversion(seqs[0], seqs[1])
        print(1.0*transition/transversion)
