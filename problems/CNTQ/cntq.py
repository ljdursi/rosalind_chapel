#!/usr/bin/env python3
from __future__ import print_function
import argparse
import sys
import re
from collections import namedtuple


Token = namedtuple('Token', ['ttype', 'name'])


class TokenType(object):
    openp = 1
    closep = 2
    done = 3
    comma = 4
    name = 5


def tokenize(newick):
    tokens = re.compile(r'(;|,|\(|\)[^\(\),;]*|[^\(\),;]+)')
    tokens = [token for token in tokens.split(newick) if token]
    token_tuples = []
    for _, token in enumerate(tokens):
        if token == "(":
            newtoken = Token(TokenType.openp, None)
        elif token[0] == ")":
            name = token[1:]
            newtoken = Token(TokenType.closep, name)
        elif token[0] == ",":
            continue
        elif token[0] == ";":
            break
        else:
            newtoken = Token(TokenType.name, token)
        token_tuples.append(newtoken)
    return token_tuples


class Tree(object):
    def __init__(self, tokens):
        self._children = []
        self._namedprogeny = []
        self.name = tokens[-1].name
        if not self.name:
            self.name = ""
        if self.name != "":
            self._namedprogeny = [self.name]

        sublists = Tree.partition_list(tokens[1:-1])
        for sublist in sublists:
            subtree = Tree(sublist)
            self._children.append(subtree)
            self._namedprogeny += subtree.namedprogeny
            if subtree.name != "":
                self._namedprogeny += [subtree.name]

    @property
    def nchildren(self):
        return len(self._children)

    @property
    def children(self):
        return self._children

    @property
    def namedprogeny(self):
        return self._namedprogeny

    def edgesets(self, child, allleafs):
        otherchildren = [c for c in self._children if c != child]
        if child.nchildren == 0:
            childsets = [set(), set()]
        else:
            childsets = [set(c.namedprogeny) for c in child.children]

        othersets = [set(c.namedprogeny) for c in otherchildren]
        if len(othersets) == 1:
            othersets.append(set(allleafs) - othersets[0] - childsets[0] - childsets[1])

        return childsets, othersets

    def edgewalk(self):
        for child in self._children:
            yield self, child
            yield from child.edgewalk()

    def nquartets(self, modulus):
        nq = 0
        leafs = set(self._namedprogeny)

        for edgep, edgec in self.edgewalk():
            childsets, othersets = edgep.edgesets(edgec, leafs)
            childtot = len(childsets[0]) + len(childsets[1])
            othertot = len(othersets[0]) + len(othersets[1])

            nlocalquartets = len(childsets[0])*len(childsets[1])*(othertot*(othertot-1))//2
            nlocalquartets += len(othersets[0])*len(othersets[1])*(childtot*(childtot-1))//2

            nq += nlocalquartets 
            
        nq = nq // 2
        return nq % modulus

    @staticmethod
    def partition_list(tokens):
        dividers = []
        level = 0
        start = None
        for i, token in enumerate(tokens):
            if token.ttype == TokenType.openp:
                level += 1
                if level == 1:
                    start = i
            elif token.ttype == TokenType.closep:
                level -= 1
                if level == 0:
                    dividers.append((start, i+1))
                    start = None
            elif token.ttype == TokenType.name and level == 0:
                dividers.append((i, i+1))
        return [tokens[s:e] for s, e in dividers]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('--modulo', '-m', type=int, default=1000000)
    args = parser.parse_args()

    with args.infile as f:
        lines = [line.strip() for line in f.readlines()]
        tree_tokens = tokenize(lines[1])
        phylogeny = Tree(tree_tokens)

        nquart = phylogeny.nquartets(args.modulo)
        print(nquart)
