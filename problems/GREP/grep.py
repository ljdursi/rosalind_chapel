#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict


def debruijn_graph(kmers):
    prefixes = defaultdict(set)

    for idx, sequence in enumerate(kmers):
        prefix = sequence[:-1]
        prefixes[prefix].add(idx)

    edges = defaultdict(set)
    for idx, sequence in enumerate(kmers):
        suffix = sequence[1:]
        for node in prefixes[suffix]:
            edges[idx].add(node)

    return edges


def dfs_paths(graph, start, sequences):
    stack = [[start]]
    while stack:
        path = stack.pop()
        vertex = path[-1]
        neighbors = graph[vertex] - set(path)
        for next_node in neighbors:
            stack.append(path + [next_node])

        if len(path) == len(sequences):
            result = "".join([sequences[idx][0] for idx in path])
            yield result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        reads = []
        for line in f:
            reads.append(line.strip())

        graph = debruijn_graph(reads)

        paths = set(dfs_paths(graph, 0, reads))
        for path in paths:
            print(path)

