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


def delta_base_counts(sequence, base1, base2):
    n = len(sequence)
    deltas = numpy.zeros(n, dtype=numpy.int)

    def inc(char, c1, c2):
        if char == c1:
            return +1
        elif char == c2:
            return -1
        else:
            return 0

    deltas[0] = inc(sequence[0], base1, base2)
    for i in range(1, n):
        base = sequence[i]
        deltas[i] = deltas[i-1] + inc(base, base1, base2)

    return deltas


def reverse_complement(sequence):
    matchings = {"A": "T", "C": "G", "G": "C", "T": "A"}
    rc = ""
    for base in sequence[::-1]:
        rc += matchings[base]

    return rc


def normalize(sequence):
    rc = reverse_complement(sequence)
    if rc < sequence:
        return sequence
    return rc


def known_good(sequences):
    normed = [normalize(seq) for seq in sequences]
    normed.sort()

    known_good = set()
    for i in range(1, len(normed)):
        if normed[i-1] == normed[i]:
            known_good.add(normed[i])
            known_good.add(reverse_complement(normed[i]))

    return known_good


def find_neighbour(bad, good_seqs):
    for seq in good_seqs:
        edit_dist = 0
        for badbase, goodbase in zip(bad, seq):
            if badbase != goodbase:
                edit_dist += 1
            if edit_dist > 1:
                break
        if edit_dist == 1:
            return seq
    return ""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-m', '--modulus', default=1000000, type=int)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(args.infile):
            sequences.append(sequence)

        known_good_seqs = known_good(sequences)
        for seq in sequences:
            if seq not in known_good_seqs:
                print(seq + "->" + find_neighbour(seq, known_good_seqs))
