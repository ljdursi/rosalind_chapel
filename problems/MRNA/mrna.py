#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict

def readtables(filename="codons.txt"):
    aa_counts = defaultdict(int)
    nstops = 0

    with open(filename, "r") as f:
        for line in f:
            items = line.strip().split()
            rna = items[0]
            amino = items[1]
            if amino == "Stop":
                nstops += 1
            else:
                aa_counts[amino] += 1
    return aa_counts, nstops

def ncombinations(protein, ncodons, nstops, modulus=1000000):
    ncombinations = nstops
    for aa in protein:
        ncombinations *= ncodons[aa]
        ncombinations = ncombinations % modulus
    return ncombinations

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    ncodons, nstops = readtables()
    with args.infile as f:
        lines = f.readlines()
        protein = "".join([l.strip() for l in lines])

        result = ncombinations(protein, ncodons, nstops)
        print(result)
