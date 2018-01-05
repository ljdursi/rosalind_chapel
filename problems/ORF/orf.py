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
    start = set()
    start.add('AUG')
    return codons, start, stop


def translate(rnaseq, codons, start, stop):
    ncodons = len(rnaseq)//3
    starts = []
    stops = []

    # translate each codon and find starts/stops
    aminos = []
    for i in range(ncodons):
        codon = rnaseq[3*i:3*i+3]
        if codon in stop:
            stops.append(i)
            aminos.append('@')
        else:
            aminos.append(codons[codon])
            if codon in start:
                starts.append(i)
    aminos = ''.join(aminos)

    # return proteins in ORFs in this sequence
    proteins = []
    laststop = -1
    for stop in stops:
        inframe_starts = [s for s in starts if s > laststop and s < stop]
        for start in inframe_starts:
            proteins.append(aminos[start:stop])
        laststop = stop

    return proteins


def transcribe(line):
    result = line.replace('T', 'U')
    return result


def reversecomplement(line):
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rc = ''.join([complements[c] for c in reversed(line)])
    return rc


def allproteins(dna, codons, start, stop):
    proteins = []

    # forward strand
    rna = transcribe(dna)
    for i in range(3):
        proteins = proteins + translate(rna[i:], codons, start, stop)

    # reverse strand
    rna = transcribe(reversecomplement(dna))
    for i in range(3):
        proteins = proteins + translate(rna[i:], codons, start, stop)

    return set(proteins)


def fasta_iterator(infile):
    started = False
    label = ""
    sequence = ""

    for line in infile:
        line = line.strip()
        if line[0] == '>':
            if started:
                yield label, sequence
            started = True
            label = line[1:]
            sequence = ""
        else:
            sequence += line

    if started:
        yield label, sequence


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    codons, start, stop = readtables()
    with args.infile as f:
        seq = ""
        for label, sequence in fasta_iterator(f):
            seq += sequence

        proteins = allproteins(seq, codons, start, stop)
        for protein in proteins:
            print(protein)
