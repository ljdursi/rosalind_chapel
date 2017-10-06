#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import urllib2
import re


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


def sequence_from_uniprot(uniprot_id):
    uniprot_url = "http://www.uniprot.org/uniprot/" + uniprot_id + ".fasta"
    response = urllib2.urlopen(uniprot_url)

    sequences = []
    for (seqlabel, sequence) in fasta_iterator(response):
        sequences.append(sequence)

    return sequences[0]


def motif_locations(sequence, motif_re):
    startidx = 0
    posns = []
    match = motif_re.search(sequence[startidx:])
    while match is not None:
        idx = startidx + match.start()
        posns.append(idx)
        startidx = idx + 1
        match = motif_re.search(sequence[startidx:])

    return posns


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--motif', default="N[^P][ST][^P]")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()
    motif = re.compile(args.motif)

    with args.infile as f:
        uniprot_ids = [line.strip() for line in f.readlines()]
        for uniprot_id in uniprot_ids:
            seq = sequence_from_uniprot(uniprot_id)
            posns = motif_locations(seq, motif)
            if not len(posns) == 0:
                print(uniprot_id)
                print(" ".join([str(posn+1) for posn in posns]))
