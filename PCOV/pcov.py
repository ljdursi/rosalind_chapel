#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict


def debruijn_graph(kmers):
    prefixes = defaultdict(set)

    for sequence in kmers:
        prefix = sequence[:-1]
        prefixes[prefix].add(sequence)

    edges = defaultdict(set)
    for sequence in kmers:
        suffix = sequence[1:]
        for node in prefixes[suffix]:
            edges[sequence].add(node)

    return edges


def reconstruct_circular_chromosome(graph):
    nodes = graph.keys()
    nodes.sort()
    current = nodes[0]
    k = len(current)
    chromosome = current

    for _ in nodes:
        current = next(iter(graph[current]))
        if not current:
            break
        chromosome += current[-1]

    return chromosome[:-k]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        reads = []
        for line in f:
            reads.append(line.strip())
        k = len(reads[0]) - 1

        g = debruijn_graph(reads)
        chrom = reconstruct_circular_chromosome(g)
        print(chrom)
