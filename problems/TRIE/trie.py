#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import namedtuple


TrieNode = namedtuple('TrieNode', ['num', 'edgedict'])


class Trie(object):
    def __init__(self):
        self.nnodes = 1
        self.root = TrieNode(1, {})

    def insert(self, sequence):
        node = self.root
        for base in sequence:
            if base in node.edgedict:
                node = node.edgedict[base]
            else:
                self.nnodes += 1
                newnode = TrieNode(self.nnodes, {})
                node.edgedict[base] = newnode
                node = newnode

    def edgelist(self):
        def dfs_edgelist(node):
            for base in node.edgedict:
                print(node.num, node.edgedict[base].num, base)
                dfs_edgelist(node.edgedict[base])
        dfs_edgelist(self.root)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-k', '--overlap', type=int, default=3)
    args = parser.parse_args()

    with args.infile as f:
        trie = Trie()
        for sequence in f.readlines():
            trie.insert(sequence.strip())

        trie.edgelist()
