#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict


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


def overlap_graph(sequences, mink):
    suffixes = defaultdict(list)
    nins = defaultdict(int)
    nsequences = len(sequences)

    for sequence, seqno in zip(sequences, range(nsequences)):
        for k in range(mink, len(sequence)):
            suffix = sequence[-k:]
            suffixes[suffix].append(seqno)

    overlaps = defaultdict(list)
    for sequence, seqno in zip(sequences, range(nsequences)):
        for k in range(mink, len(sequence)):
            prefix = sequence[:k]
            locations = suffixes[prefix]
            for idx in locations:
                if not idx == seqno:
                    overlaps[idx].append((seqno, k))
                    nins[seqno] += 1

    print(nins)
    startpos = [idx for idx in range(nsequences) if nins[idx] == 0]
    return overlaps, startpos


def reconstruct_chromosome(graph, startpos, seqs):
    # we are assured there's a unique reconstruction
    current = startpos[0]
    chromosome = seqs[current]

    while not len(graph[current]) == 0:
        current, overlap = graph[current][0]
        chromosome += seqs[current][overlap:]

    return chromosome


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        reads = []
        for label, sequence in fasta_iterator(f):
            reads.append(sequence)
        l = len(reads[0])
        k = (l+1)//2

        g, starts = overlap_graph(reads, k)
        chrom = reconstruct_chromosome(g, starts, reads)
        print(chrom)
