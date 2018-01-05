#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import re


class Tokens(object):
    openp = 1
    closep = 2
    close_w_parent = 3
    done = 4
    comma = 5
    label = 6


def tokenize(newick):
    tokens = re.compile('(;|,|\(|\)[^\(\),;]+|\)|[^\(\),;]+)')
    tokens = [token for token in tokens.split(newick) if token]
    tokenpairs = []
    last = len(tokens)-1
    for i, token in enumerate(tokens):
        if token == "(":
            pair = (Tokens.openp, token)
        elif token == ")":
            pair = (Tokens.closep, token)
        elif token[0] == ")":
            pair = (Tokens.close_w_parent, token[1:])
        elif token[0] == ",":
            pair = (Tokens.comma, token)
        elif token[0] == ";":
            pair = (Tokens.done, token)
            last = i
        else:
            pair = (Tokens.label, token)
        tokenpairs.append(pair)
    return tokenpairs[:last+1]


def find_dist(newick, node1, node2):
    tokens = tokenize(newick)
    loc1 = [i for i, tok in enumerate(tokens) if tok[1] == node1][0]
    loc2 = [i for i, tok in enumerate(tokens) if tok[1] == node2][0]

    start, finish = loc1, loc2
    if loc2 < loc1:
        start, finish = loc2, loc1
        node1, node2 = node2, node1

    state = [None]
    for token_type, token_val in tokens[start+1:finish+1]:
        if token_type in [Tokens.closep, Tokens.close_w_parent]:
            if state[-1] == Tokens.openp:
                state.pop()
            else:
                state.append(token_type)
        elif token_type == Tokens.openp:
            state.append(token_type)

    state = state[1:]
    if token_type == Tokens.close_w_parent:
        return len(state)
    else:
        return len(state) + 2


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
