#!/usr/bin/env python
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
    tokens = re.compile('(;|,|\(|\)[^\(\),;]*|[^\(\),;]+)')
    tokens = [token for token in tokens.split(newick) if token]
    token_tuples = []
    for token in tokens:
        if token == "(":
            tt = Token(TokenType.openp, None)
        elif token[0] == ")":
            tt = Token(TokenType.closep, token[1:])
        elif token[0] == ",":
            continue
        elif token[0] == ";":
            break
        else:
            tt = Token(TokenType.name, token)
        token_tuples.append(tt)
    return token_tuples


class Tree(object):
    def __init__(self, tokens):
        self._children = []
        self._namedprogeny = []
        self.name = tokens[-1].name
        if not self.name:
            self.name = ""

        sublists = Tree.partition_list(tokens[1:-1])
        for sublist in sublists:
            subtree = Tree(sublist)
            self._children.append(subtree)
            self._namedprogeny += subtree.namedprogeny
            if subtree.name != "":
                self._namedprogeny += [subtree.name]

    @property
    def namedchildren(self):
        return [kid.name for kid in self._children if kid.name != ""]

    @property
    def namedprogeny(self):
        return self._namedprogeny

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

    def __iter__(self):
        for child in self._children:
            for subtree in child:
                yield subtree
        yield self

    def splits(self, allnames):
        n = len(allnames)
        allnames_set = set(allnames)
        split_set = set([])

        for subtree in self:
            subnames = set(subtree.namedprogeny)
            if len(subnames) in [0, 1, n-1, n]:
                continue
            if len(subnames) > n-len(subnames):
                subnames = allnames_set - subnames

            split_set.add(tuple(subnames))
        return split_set


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = [line.strip() for line in f.readlines()]
        allnames = [item.strip() for item in lines[0].split()]
        phylogeny1 = Tree(tokenize(lines[1]))
        phylogeny2 = Tree(tokenize(lines[2]))

        splitset1 = phylogeny1.splits(allnames)
        splitset2 = phylogeny2.splits(allnames)

        print(len(splitset1 ^ splitset2))
