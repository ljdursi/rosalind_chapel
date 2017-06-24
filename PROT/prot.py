#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys

def readtables(filename="codons.txt"):
    codons = {}
    stop = set()

    with open(filename, "r") as f:
        for line in f:
            items = line.strip().split()
            rna = items[0]
            amino = items[1]
            if amino == "Stop":
                stop.add(rna)
            else:
                codons[rna] = amino
    return codons, stop

def translate(rnaseq, codons, stop):
    ncodons = len(rnaseq)//3
    protein = ""
    for i in range(ncodons):
        codon = rnaseq[3*i:3*i+3]
        if codon in stop:
            return protein
        else:
            protein = protein + codons[codon]
    return protein

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    codons, stop = readtables()
    with args.infile as f:
        lines = f.readlines()
        rnaseq = "".join([l.strip() for l in lines])

        result = translate(rnaseq, codons, stop)
        print(result)
