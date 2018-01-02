#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys


def char_table_from_strings(lines):
    nlines = len(lines)
    nnucleotides = len(lines[0])

    characters = []
    for i in range(nnucleotides):
        character = "1"
        one_base = lines[0][i]
        for j in range(1, nlines):
            base = lines[j][i]
            if base == one_base:
                character += "1"
            else:
                character += "0"

        characters.append(character)
    return characters


def trivial(character):
    ones = character.count("1")
    zeros = character.count("0")
    if ones in [0, 1] or zeros in [0, 1]:
        return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]

        ctbl = char_table_from_strings(lines)
        for character in ctbl:
            if not trivial(character):
                print(character)
