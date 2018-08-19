#!/usr/bin/env python
#pylint: disable=W1401
#pylint: disable=C0103
"""
Alignment-Based Phylogeny
http://rosalind.info/problems/alph/
"""
from __future__ import print_function
import argparse
import sys
import re
from collections import namedtuple
import numpy

_BASE_TO_IDX = {'A':0, 'C':1, 'G':2, 'T':3, '-':4}
_NBASES = len(_BASE_TO_IDX)
_IDX_TO_BASE = {val:key for key, val in _BASE_TO_IDX.items()}

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
    for _, token in enumerate(tokens):
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


def hamming(sequence1, sequence2):
    """
    Hamming distance between two strings
    """
    score = 0
    assert len(sequence1) == len(sequence2)
    for (base1, base2) in zip(sequence1, sequence2):
        if base1 != base2:
            score += 1
    return score


class Tree(object):
    def __init__(self, tokens, sequences, k):
        self._children = []
        self.name = tokens[-1].name
        self.sequence = sequences[self.name] if self.name in sequences else None

        sublists = Tree.partition_list(tokens[1:-1])
        for sublist in sublists:
            subtree = Tree(sublist, sequences, k)
            self._children.append(subtree)

        if self.sequence:
            self.seqmatrix = self._sequence_to_matrix()
        else:
            self.seqmatrix = sum([self._project_matrix(child.seqmatrix) for child in self._children])

    def _sequence_to_matrix(self):
        matrix = numpy.zeros((len(self.sequence), _NBASES), dtype=numpy.int)
        matrix += 1
        for idx, base in enumerate(self.sequence):
            matrix[idx, _BASE_TO_IDX[base]] = 0
        return matrix

    def _project_matrix(self, mtx):
        k, _ = mtx.shape
        penalty = 1 - numpy.eye(_NBASES, dtype=numpy.int)
        newmtx = numpy.zeros_like(mtx)
        for i in range(k):
            row = mtx[i, :] 
            choices = numpy.repeat(row, _NBASES).reshape((_NBASES, _NBASES)) + penalty
            newmtx[i, :] = numpy.amin(choices, axis=0)
        return newmtx

    def _seqmatrix_to_sequence(self):
        if self.sequence:
            return
        self.sequence = ""
        k, _ = self.seqmatrix.shape
        for i in range(k):
            self.sequence += _IDX_TO_BASE[numpy.argmin(self.seqmatrix[i, :])]

    def set_sequences(self):
        self._seqmatrix_to_sequence()
        matrix = self._sequence_to_matrix()
        for child in self._children:
            child.seqmatrix += matrix
            child.set_sequences()

    def calc_hamming(self):
        tree_hamming = 0
        for child in self._children:
            tree_hamming += child.calc_hamming()
            tree_hamming += hamming(self.sequence, child.sequence)

        return tree_hamming

    def sequences(self):
        sequences_dict = {self.name: self.sequence}
        for child in self._children:
            child_dict = child.sequences()
            sequences_dict.update(child_dict)
        return sequences_dict

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


def fasta_to_dict(lines):
    """
    Simple fasta parser
    """
    started = False
    label = ""
    sequence = ""
    sequence_dict = {}

    for line in lines:
        line = line.strip()
        if line[0] == '>':
            if started:
                sequence_dict[label] = sequence
            started = True
            label = line[1:]
            sequence = ""
        else:
            sequence += line

    if started:
        sequence_dict[label] = sequence

    return sequence_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = [line.strip() for line in f.readlines()]
        tokens = tokenize(lines[0])
        given_sequences = fasta_to_dict(lines[1:])

        k = None
        for _, seq in given_sequences.items():
            if not k:
                k = len(seq)
            else:
                assert len(seq) == k

        phylogeny = Tree(tokens, given_sequences, k)
        phylogeny.set_sequences()

        tree_sequences = phylogeny.sequences()
        score = phylogeny.calc_hamming()

        print(score)
        for taxon, seq in tree_sequences.items():
            if not taxon in given_sequences:
                print(">"+taxon)
                print(seq)

if __name__ == "__main__":
    main()
