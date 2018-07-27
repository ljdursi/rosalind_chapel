#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict


def debruijn_graph(kmers):
    prefixes = defaultdict(set)

    for idx, sequence in enumerate(kmers):
        prefix = sequence[:-1]
        prefixes[prefix].add((idx, sequence))

    edges = defaultdict(set)
    for idx, sequence in enumerate(kmers):
        suffix = sequence[1:]
        for node in prefixes[suffix]:
            edges[(idx, sequence)].add(node)

    return edges


def dfs_paths(graph, start):
    stack = [(start, [start])]
    while stack:
        (vertex, path) = stack.pop()
        for next_node in graph[vertex] - set(path):
            stack.append((next_node, path + [next_node]))

        yield "".join([seq[0] for idx, seq in path])


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
        paths = set(dfs_paths(graph, (0, reads[0])))
        fullpaths = [path for path in paths if len(path) == len(reads)]
        for path in fullpaths:
            print(path)

