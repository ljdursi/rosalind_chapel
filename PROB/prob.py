#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
from collections import defaultdict
from math import log10


def stringprobs(seq, gcprobs):
    nbases = defaultdict(int)
    for base in seq:
        nbases[base] += 1

    logprobs = []
    for gcprob in gcprobs:
        atprob = 1. - gcprob
        logaprob = log10(atprob/2.)
        loggprob = log10(gcprob/2.)

        logprob = nbases['G']*loggprob + nbases['C']*loggprob
        logprob += nbases['A']*logaprob + nbases['T']*logaprob

        logprobs.append(logprob)

    return logprobs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]

        sequence = lines[0]
        gcprobs = [float(item) for item in lines[1].split()]
        sequenceprobs = stringprobs(sequence, gcprobs)
        print(" ".join([str(prob) for prob in sequenceprobs]))
