#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys

def readtable(filename="weights.txt"):
    weights = {}

    with open(filename, "r") as f:
        for line in f:
            items = line.strip().split()
            amino = items[0]
            weight = float(items[1])
            weights[amino] = weight
    return weights

def totalweight(protein, weights):
    weight = 0.0
    for amino in protein:
        weight += weights[amino]
    return weight

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    weights = readtable()
    with args.infile as f:
        lines = f.readlines()
        protein = "".join([l.strip() for l in lines])

        result = totalweight(protein, weights)
        print(result)
