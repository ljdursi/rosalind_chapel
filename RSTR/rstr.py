#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
from collections import defaultdict
import numpy


def string_logprob(seq, gcprob):
    nbases = defaultdict(int)
    for base in seq:
        nbases[base] += 1

    atprob = 1. - gcprob
    logaprob = numpy.log(atprob/2.)
    loggprob = numpy.log(gcprob/2.)

    logprob = nbases['G']*loggprob + nbases['C']*loggprob
    logprob += nbases['A']*logaprob + nbases['T']*logaprob

    return logprob


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        items = lines[0].split()
        n = int(items[0])
        x = float(items[1])

        sequence = lines[1]
        sequence_non_prob = -numpy.expm1(string_logprob(sequence, x))
        p = 1. - sequence_non_prob**n
        print(p)
