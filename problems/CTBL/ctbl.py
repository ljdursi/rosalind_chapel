#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import re
from collections import namedtuple


Token = namedtuple('Token', ['ttype', 'name', 'weight'])


class TokenType(object):
    openp = 1
    closep = 2
    done = 3
    comma = 4
    name = 5


def parseweight(name):
    if len(name) == 0:
        return None, 1
    items = name.split(':')
    if len(items) == 1:
        return items[0], 1
    else:
        return items[0], int(items[1])


def tokenize(newick):
    tokens = re.compile('(;|,|\(|\)[^\(\),;]*|[^\(\),;]+)')
    tokens = [token for token in tokens.split(newick) if token]
    token_tuples = []
    for i, token in enumerate(tokens):
        if token == "(":
            tt = Token(TokenType.openp, None, None)
        elif token[0] == ")":
            name, weight = parseweight(token[1:])
            tt = Token(TokenType.closep, name, weight)
        elif token[0] == ",":
            continue
        elif token[0] == ";":
            break
        else:
            name, weight = parseweight(token)
            tt = Token(TokenType.name, name, weight)
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


def character_entry(allnames, subnames):
    char_in, char_out = '1', '0'
    if subnames[0] == allnames[0]:
        char_in, char_out = '0', '1'

    entry = ""
    for name in allnames:
        if name in subnames:
            entry += char_in
        else:
            entry += char_out

    return entry


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            tokens = tokenize(line)
            phylogeny = Tree(tokens)

            allnames = sorted(phylogeny.namedprogeny)
            nnames = len(allnames)
            for i, subtree in enumerate(phylogeny):
                subnames = sorted(subtree.namedprogeny)
                if len(subnames) in [0, 1, nnames-1, nnames]:
                    continue
                print(character_entry(allnames, subnames))
