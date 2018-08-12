#!/usr/bin/env python
"""
Enumerating Unrooted Binary Trees
http://rosalind.info/problems/EUBT/
"""
from __future__ import print_function
import argparse
import sys
from collections import namedtuple
from copy import deepcopy

Edge = namedtuple("Edge", ["node1idx", "node2idx"])

class Tree(object):
    """
    Unrooted binary tree
    """

    class Node(object):
        """
        Node in an unrooted binary tree
        """
        def __init__(self, label, idx):
            self.label = label
            self.edgeidxs = []
            self.idx = idx

        @property
        def nedges(self):
            return len(self.edgeidxs)

        def add_edge(self, edgeidx):
            self.edgeidxs.append(edgeidx)

        def remove_edge(self, edgeidx):
            if edgeidx not in self.edgeidxs:
                return
            ei = self.edgeidxs.index(edgeidx)
            self.edgeidxs[ei] = self.edgeidxs[-1]
            self.edgeidxs.pop()


    def __init__(self, nodelabels, edges):
        self.nodes = []
        self.edges = []

        for lab in nodelabels:
            self.add_node(lab)

        for edge in edges:
            self.add_edge(edge[0], edge[1])

    @property
    def num_nodes(self):
        return len(self.nodes)

    @property
    def num_edges(self):
        return len(self.edges)

    def add_node(self, nodelabel):
        nodenum = self.num_nodes
        self.nodes.append(self.Node(nodelabel, nodenum))
        return nodenum

    def add_edge(self, node1idx, node2idx):
        edgenum = self.num_edges
        edge = Edge(node1idx=node1idx, node2idx=node2idx)
        self.edges.append(edge)
        self.nodes[node1idx].add_edge(edgenum)
        self.nodes[node2idx].add_edge(edgenum)
        return edgenum

    def _split_edge(self, edgeidx, nodelabel):
        newidx = self.add_node(nodelabel)
        (n1, n2) = self.edges[edgeidx]
        self.nodes[n2].remove_edge(edgeidx)
        self.edges[edgeidx] = Edge(n1, newidx)
        self.nodes[newidx].add_edge(edgeidx)

        self.add_edge(newidx, n2)
        return newidx

    def insert_node_on_edge(self, edgeidx, nodelabel):
        internal_node = self._split_edge(edgeidx, "")
        newnodeidx = self.add_node(nodelabel)
        self.add_edge(internal_node, newnodeidx)

    def node_neighbours(self, nodeidx):
        neighs = []
        node = self.nodes[nodeidx]

        for i in range(node.nedges):
            edgeidx = node.edgeidxs[i]
            neighbour = self.edges[edgeidx].node1idx
            if neighbour == nodeidx:
                neighbour = self.edges[edgeidx].node2idx
            neighs.append(neighbour)

        return neighs

    def _newick_dfs(self, nodeidx, visited):
        subtreestr = ""

        if visited[nodeidx] or nodeidx < 0 or nodeidx > self.num_nodes:
            return subtreestr

        visited[nodeidx] = True
        node = self.nodes[nodeidx]
        open_neighbours = [neigh for neigh in self.node_neighbours(nodeidx) if not visited[neigh]]

        if len(open_neighbours) == 0:
            return node.label

        subtreestr = "("
        for neigh in open_neighbours[:-1]:
            subtreestr = subtreestr + self._newick_dfs(neigh, visited) + ","

        subtreestr = subtreestr + self._newick_dfs(open_neighbours[-1], visited) + ")"
        subtreestr = subtreestr + node.label

        return subtreestr

    def newick(self):
        if self.num_nodes == 0:
            return ""

        visited = [False]*self.num_nodes
        return self._newick_dfs(0, visited)


def enumerate_trees(base_tree, newnodes):
    if len(newnodes) == 0:
        print(base_tree.newick()+";")
        return

    for edgeidx in range(base_tree.num_edges):
        subtree = deepcopy(base_tree)
        subtree.insert_node_on_edge(edgeidx, newnodes[0])
        enumerate_trees(subtree, newnodes[1:])

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            all_names = line.split()

            tree = Tree(all_names[0:2], [(0, 1)])
            enumerate_trees(tree, all_names[2:])

