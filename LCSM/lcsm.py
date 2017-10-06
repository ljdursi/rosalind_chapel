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


def substr_iterator(sequence):
    l = len(sequence)
    for sublen in range(l, 0, -1):
        for start in range(0, l-sublen+1):
            yield sequence[start:start+sublen]


def naive_lcsm(sequences):
    sorted_seqs = sorted(sequences, key=len)
    shortest = sorted_seqs[0]
    others = sorted_seqs[1:]

    for substr in substr_iterator(shortest):
        allmatch = True
        for otherstr in others:
            if substr not in otherstr:
                allmatch = False
                break
        if allmatch:
            return substr

    return ""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(f):
            sequences.append(sequence)

        lcsm = naive_lcsm(sequences)
        print(lcsm)
