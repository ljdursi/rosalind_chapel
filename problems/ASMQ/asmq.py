#!/usr/bin/env python
"""
Assessing Assembly Quality with N50 and N75
http://rosalind.info/problems/asmq/
"""
from __future__ import print_function
import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        read_lengths = []
        for line in f:
            read_lengths.append(len(line.strip()))

        n_tot = sum(read_lengths)
        read_lengths.sort()

        tot_so_far = 0
        for rlen in read_lengths:
            if (n_tot - tot_so_far) >= 3*n_tot//4:
                n75 = rlen
            if (n_tot - tot_so_far) >= n_tot//2:
                n50 = rlen
            tot_so_far += rlen

        print(n50, n75)
