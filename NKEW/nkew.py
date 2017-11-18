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
    last = len(tokens)-1
    for i, token in enumerate(tokens):
        if token == "(":
            tt = Token(TokenType.openp, None, None)
        elif token[0] == ")":
            name, weight = parseweight(token[1:])
            tt = Token(TokenType.closep, name, weight)
        elif token[0] == ",":
            tt = Token(TokenType.comma, None, None)
        elif token[0] == ";":
            tt = Token(TokenType.done, None, None)
            last = i
        else:
            name, weight = parseweight(token)
            tt = Token(TokenType.name, name, weight)
        token_tuples.append(tt)
    return token_tuples[:last+1]


def find_dist(newick, node1, node2):
    tokens = tokenize(newick)

    loc1, loc2 = 0, 0
    for i, token in enumerate(tokens):
        if token.name == node1:
            loc1 = i
        if token.name == node2:
            loc2 = i

    start, finish = loc1, loc2
    if loc1 > loc2:
        node1, node2 == node2, node1
        start, finish = loc2, loc1

    level, minlevel = 0, 0
    path = [tokens[start]]

    for idx in range(start+1, len(tokens)):
        token = tokens[idx]
        if token.ttype == TokenType.closep:
            level -= 1
            if path[-1][0] == TokenType.openp:
                path.pop()
            else:
                path.append(token)
        elif token.ttype == TokenType.openp:
            level += 1
            path.append(token)
        elif token.name in [node1, node2]:
            path.append(token)

        # we may have to keep seeing )s to find remaining edge weights
        if level < minlevel and idx < finish:
            minlevel = level
        if idx >= finish and level <= minlevel:
            break

    pathlen = sum([t.weight for t in path if t.weight])
    if path[-1].name == node2 and level < minlevel:
        pathlen -= path[-1].weight

    return pathlen


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines]
        lines = [l for l in lines if len(l) != 0]

        dists = []
        for newick, itemline in zip(lines[::2], lines[1::2]):
            items = itemline.split()
            dists.append(find_dist(newick, items[0], items[1]))

        print(" ".join([str(d) for d in dists]))
