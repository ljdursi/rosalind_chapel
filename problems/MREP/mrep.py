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
            self.label = ""
            self.left_char = None

        def __repr__(self):
            return str(self.left_char) + " " + self.label + " " + str(len(self._edges)) 

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

        def printnodes(self, prefix=""):
            print(prefix + self.__repr__())
            for _, v in self._edges.items():
                v.child.printnodes(prefix+"  ")
    
        def set_left_chars(self):
            for _, v in self._edges.items():
                v.child.set_left_chars()
            left_chars = [v.child.left_char for _, v in self._edges.items() if v.child]
            if len(left_chars) > 0:
                if all([item == left_chars[0] for item in left_chars]):
                    self.left_char = left_chars[0]
                else:
                    self.left_char = "*"

        def maximal_repeats(self):
            repeats = []
            children = [ v.child for _, v in self._edges.items() if v.child]
            if len(children) >= 2:
                for child in children:
                    repeats += child.maximal_repeats()
                if self.left_char == '*':
                    repeats += [self.label]
            return repeats

    class Edge(object):
        def __init__(self, label, parent, child):
            self._parent = parent
            self._child = child
            self._label = label
            self._len = len(label)

        def __repr__(self):
            rep = "Edge(label=" + self._label
            rep += ", child=" + self._child.__repr__() + ")"
            return rep

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

            plabel = "" if not self._parent else self._parent.label
            newnode.label = plabel + label[:idx]
            return newnode

    def _insert(self, string, left_char="#"):
        node, edgelab, ncommon, sremaining = self._root.walk(string)

        if sremaining == "":
            return

        if edgelab is not None:
            edge = node.edge(edgelab)
            node = edge.splitAt(ncommon)

        newnode = NaiveQuadratic.Node()
        newnode.label = string
        newnode.left_char = left_char

        newedge = NaiveQuadratic.Edge(sremaining, node, newnode)
        node.add_out_edge(newedge)
        return

    def __init__(self, string):
        self._root = NaiveQuadratic.Node()
        for i in range(len(string)):
            self._insert(string[i:], None if i == 0 else string[i-1])
        self._root.set_left_chars()

    def printedges(self):
        self._root.printedges()

    def printnodes(self):
        self._root.printnodes()

    def maximal_repeats(self):
        return self._root.maximal_repeats()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-k', '--threshold', type=int, default=20)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            stree = NaiveQuadratic(line.strip()+"$")
            for rep in stree.maximal_repeats():
                if len(rep) >= args.threshold:
                    print(rep)
