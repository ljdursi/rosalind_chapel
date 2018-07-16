#!/usr/bin/env python
"""
Genome Assembly With Reads
http://rosalind.info/problems/gasm/
"""
from __future__ import print_function
import argparse
import sys
from collections import defaultdict

_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def reverse_complement(seq):
    """
    Reverse complement of a sring
    """
    rc = ""
    for base in seq[::-1]:
        rc += _COMPLEMENT[base]
    return rc


def cyclic_rc_match(seq1, seq2):
    """
    Are these strings consistent with being cyclic
    shifts of each other, one being a reverse
    complement?
    """
    if len(seq1) != len(seq2):
        return False
    return reverse_complement(seq2) in seq1+seq1


def generate_kmers(reads, k):
    """
    The list of kmers, possibly not unique, that
    can come from the given reads or their reverse
    complements
    """
    kmers = []
    n = len(reads[0])
    for read in reads:
        for i in range(n-k+1):
            kmers.append(read[i:i+k])
        rc_read = reverse_complement(read)
        for i in range(n-k+1):
            kmers.append(rc_read[i:i+k])
    return kmers
    

def debruijn_graph(kmers):
    """
    Return a dictionary of edges:
    start_vtx => [end_vtx1, end_vtx2...]
    """
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
    chromosomes = []
    k = None

    def get_next(key):
        if not key in graph:
            return None
        try:
            nextnode = graph[key].pop()
            if len(graph[key]) == 0:
                graph.pop(key)
        except:
            nextnode = None
        return nextnode

    while len(graph) > 0:
        nodes = graph.keys()
        nodes.sort()
        if len(nodes) == 0:
            break
        current = nodes[0]
        if not k:
            k = len(current)

        cur_chromosome = current
        nextnode = get_next(current)

        while nextnode:
            cur_chromosome += nextnode[-1]
            nextnode = get_next(nextnode)

        chromosomes.append(cur_chromosome[:-k])

    return chromosomes


def success(chroms):
    if len(chroms) != 2:
        return False
    return cyclic_rc_match(chroms[0], chroms[1])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        reads = []
        for line in f:
            reads.append(line.strip())

        for k in range(2, len(reads[0]))[::-1]:
            kmers = generate_kmers(reads, k)
            g = debruijn_graph(kmers)
            chromosomes = reconstruct_circular_chromosome(g)
            if success(chromosomes):
                print(chromosomes[0])
                break
