#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import itertools
from sets import Set


def partition_from_character(character):
    part_zero = []
    part_one = []
    for i, ch in enumerate(character):
        if ch == "0":
            part_zero.append(i)
        elif ch == "1":
            part_one.append(i)

    return part_zero, part_one


def quartets_from_partition(seta, setb):
    for first in itertools.combinations(seta, 2):
        for second in itertools.combinations(setb, 2):
            if first < second:
                yield (first[0], first[1], second[0], second[1])
            else:
                yield (second[0], second[1], first[0], first[1])


def quartet_string(allnames, q):
    result = "{" + allnames[q[0]] + ", " + allnames[q[1]] + "} "
    result += "{" + allnames[q[2]] + ", " + allnames[q[3]] + "}"
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    quartets = Set()
    with args.infile as f:
        lines = [line.strip() for line in f.readlines()]
        allnames = lines[0].split()
        for line in lines[1:]:
            part_zero, part_one = partition_from_character(line)
            for quartet in quartets_from_partition(part_zero, part_one):
                quartets.add(quartet)
        for quartet in quartets:
            print(quartet_string(allnames, quartet))
