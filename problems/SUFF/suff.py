#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import os.path

class NaiveQuadratic(object):
    """
    Naive method: uses quadratic space, quadratic build time
    Inserts all suffixes in decreasing order of length
    """
    class Node(object):
        def __init__(self):
            self._edges = {}

        def edge(self, key):
            return self._edges[key]

        def add_out_edge(self, edge):
            start = edge.label[0]
            if start in self._edges.keys():
                raise KeyError("Edge starting with "+start+" already in node")
            self._edges[start] = edge

        def walk(self, string):
            """
            Walks the tree matching a string as far as it can.
            Returns the last node found on the path, the first character
            of the edge (if such an edge exists), the index into
            the edge and the remaining string (respectively) where the match
            stops.
            """
            if string == "" or string[0] not in self._edges.keys():
                return self, None, None, string
            edge = self._edges[string[0]]

            common = os.path.commonprefix([string, edge.label])
            ncommon = len(common)

            if ncommon == edge.length:
                return edge.child.walk(string[ncommon:])

            return self, string[0], ncommon, string[ncommon:]

        def printedges(self):
            for k, v in self._edges.items():
                print(v.label)
                if v.child:
                    v.child.printedges()

    class Edge(object):
        def __init__(self, label, parent, child):
            self._parent = parent
            self._child = child
            self._label = label
            self._len = len(label)

        @property
        def label(self):
            return self._label

        @property
        def length(self):
            return self._len

        @property
        def child(self):
            return self._child

        def splitAt(self, idx):
            if idx >= self._len:
                raise IndexError("Split index out of range")

            label = self._label
            newnode = NaiveQuadratic.Node()
            newedge = NaiveQuadratic.Edge(label[idx:], newnode, self._child)
            newnode.add_out_edge(newedge)
            self._label = label[:idx]
            self._len = len(self._label)
            self._child = newnode
            return newnode

    def _insert(self, string):
        node, edgelab, ncommon, sremaining = self._root.walk(string)

        if sremaining == "":
            return

        if edgelab is not None:
            edge = node.edge(edgelab)
            node = edge.splitAt(ncommon)

        newedge = NaiveQuadratic.Edge(sremaining, node, None)
        node.add_out_edge(newedge)
        return

    def __init__(self, string):
        self._root = NaiveQuadratic.Node()
        for i in range(len(string)):
            self._insert(string[i:])

    def printedges(self):
        self._root.printedges()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            stree = NaiveQuadratic(line.strip())
            stree.printedges()
