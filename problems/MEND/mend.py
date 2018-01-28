#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import re
import numpy
from collections import namedtuple


def allele_count_probabilities(parent1, parent2):
    # calculate transition matrices
    T0 = numpy.zeros((3, 3))
    T2 = numpy.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            T0[i, j] = (2 - i)*(2 - j)/4.0
            T2[i, j] = i*j/4.0

    allele_probs = numpy.zeros(3)
    allele_probs[0] = numpy.dot(parent1, numpy.dot(T0, parent2))
    allele_probs[2] = numpy.dot(parent1, numpy.dot(T2, parent2))
    allele_probs[1] = 1. - allele_probs[0] - allele_probs[2]
    return allele_probs


def prob_vector_from_genotype(genotype):
    prob_dict = {'AA': (1., 0., 0.), 'Aa': (0., 1., 0.),
                 'aA': (0., 1., 0.), 'aa': (0., 0., 1.)}
    return numpy.array(prob_dict[genotype])


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
    for i, token in enumerate(tokens):
        if token == "(":
            tt = Token(TokenType.openp, None)
        elif token[0] == ")":
            name = token[1:]
            tt = Token(TokenType.closep, name)
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
        self._gt_probs = None
        if not self.name:
            self.name = ""
        else:
            self._gt_probs = prob_vector_from_genotype(self.name)

        sublists = Tree.partition_list(tokens[1:-1])
        for sublist in sublists:
            subtree = Tree(sublist)
            self._children.append(subtree)
            if subtree.name != "":
                self._namedprogeny += [subtree.name]

        if self._gt_probs is None:
            if len(self._children) == 2:
                gtps = [c.gt_probability for c in self._children]
                self._gt_probs = allele_count_probabilities(gtps[0], gtps[1])

    @property
    def gt_probability(self):
        return self._gt_probs

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


def n_choose_k_mod_m(n, k, modulus):
    if k < 0 or k > n:
        return 0

    if n == k or k == 0:
        return 1

    num, den = 1, 1
    for i in range(1, min(k, n-k) + 1):
        num *= n
        den *= i
        n -= 1
    return (num // den) % modulus


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = [line.strip() for line in f.readlines()]
        tokens = tokenize(lines[0])
        phylogeny = Tree(tokens)

        print(" ".join([str(p) for p in phylogeny.gt_probability]))
