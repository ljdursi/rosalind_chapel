#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
from collections import defaultdict


class Tree(object):
    @classmethod
    def partition_list(cls, lines, children):
        sublists = defaultdict(list)
        subtrees = {}

        for i, child in enumerate(children):
            subtrees[child] = i

        for line in lines:
            items = line.split()
            parent_name, child_name = items[0], items[1]
            if parent_name in subtrees:
                subtree_num = subtrees[parent_name]
                sublists[subtree_num].append(line)
                subtrees[child_name] = subtree_num

        return [sublists[i] for i in range(len(children))]

    def __init__(self, lines, fullstring, name,
                 start=0, length=0, totlength=0):
        self.l = len(fullstring)
        self.children = []
        self.name = name
        self.start = start
        self.length = length
        self.totlength = totlength

        self.nendings = 0
        if start + length >= self.l:
            self.nendings = 1

        if len(lines) == 0:
            return

        children = []
        i = 0
        for line in lines:
            items = line.split()
            if items[0] == name:
                children.append((items[1], int(items[2]), int(items[3])))
                i += 1
            else:
                break
        lines = lines[i:]

        sublists = Tree.partition_list(lines, [item[0] for item in children])
        for (sublist, child) in zip(sublists, children):
            childname, childstart, childlen = child
            childtotlen = childlen + totlength
            child = Tree(sublist, fullstring, childname, childstart,
                         childlen, childtotlen)
            self.nendings += child.nendings
            self.children.append(child)

    def __iter__(self):
        for child in self.children:
            for subtree in child:
                yield subtree
        yield self


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = [line.strip() for line in f.readlines()]
        superseq = lines[0]
        minrepeats = int(lines[1])
        treelines = lines[2:]

        rootname = treelines[0].split()[0]
        suffixtree = Tree(treelines, superseq, rootname)

        maxlen, maxstart, maxname = 0, 0, ""
        for subtree in suffixtree:
            if subtree.nendings >= minrepeats:
                if subtree.totlength > maxlen:
                    maxlen, maxname = subtree.totlength, subtree.name
                    maxstart = subtree.start - subtree.totlength
                    maxstart += subtree.length
        print(superseq[maxstart-1:maxstart-1+maxlen])
