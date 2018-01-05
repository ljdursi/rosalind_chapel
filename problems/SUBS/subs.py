#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        seq = lines[0].strip()
        nseq = len(seq)
        query = lines[1].strip()
        nquery = len(query)

        locations = []
        for i in range(nseq-nquery+1):
            if seq[i:i+nquery] == query:
                locations.append(i+1)
        print(" ".join([str(l) for l in locations]))
