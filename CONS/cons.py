#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import numpy

base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
idx_to_base = {}
for key in base_to_idx:
    idx_to_base[base_to_idx[key]] = key

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

def profile(sequences):
    l = len(sequences[0])
    profile = numpy.zeros((4, l), numpy.uint32)
    for sequence in sequences:
        for i in range(l):
            base = sequence[i]
            profile[base_to_idx[base], i] += 1
    return profile

def consensus_from_profile(profile):
    indices = numpy.argmax(profile, axis=0)
    consensus = "".join([idx_to_base[index] for index in indices])
    return consensus

def output_profile(profile):
    nrows, ncols = profile.shape
    for r in range(nrows):
        print(idx_to_base[r]+ ": ", end='')
        print(" ".join([str(count) for count in profile[r,:]]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        sequences = []
        for label, sequence in fasta_iterator(f):
            sequences.append(sequence)
        counts = profile(sequences)
        consensus = consensus_from_profile(counts)
        print(consensus)
        output_profile(counts)
