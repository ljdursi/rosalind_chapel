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

    return overlaps


def dfs_traversals(graph, start, curvisited):
    visited = curvisited.copy()
    visited.add(start)
    neighbors = [(i, k) for i, k in graph[start] if i not in visited]
    if len(neighbors) == 0:
        return [[(start, None)]]
    else:
        traversals = []
        for neighbor, overlap in neighbors:
            subtraversals = dfs_traversals(graph, neighbor, visited)
            traversals += [[(start, overlap)] + sub for sub in subtraversals]
    return traversals


def possible_chromosomes(graph, seqs):
    nseqs = len(seqs)
    nins = [0]*nseqs
    for k, values in graph.items():
        for idx, k in values:
            nins[idx] += 1

    possible_starts = [i for i, j in zip(range(nseqs), nins) if j == 0]

    full_traversals = []
    for start in possible_starts:
        visited = set()
        traversals = dfs_traversals(graph, start, visited)
        for traversal in traversals:
            if len(traversal) == nseqs:
                full_traversals.append(traversal)

    chromosomes = []
    for traversal in full_traversals:
        idxs = [i for i, j in traversal]
        overlaps = [0] + [j for i, j in traversal[:-1]]

        chromosome = "".join([seqs[idx][overlap:] for idx, overlap in zip(idxs, overlaps)])
        chromosomes.append(chromosome)

    return chromosomes


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

        g = overlap_graph(reads, k)
        chroms = possible_chromosomes(g, reads)

        minchrom = min(chroms, key=len)
        print(minchrom)
