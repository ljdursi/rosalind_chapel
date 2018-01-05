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
    logaprob = numpy.log(atprob/2.+1.e-9)
    loggprob = numpy.log(gcprob/2.+1.e-9)

    logprob = nbases['G']*loggprob + nbases['C']*loggprob
    logprob += nbases['A']*logaprob + nbases['T']*logaprob

    return logprob


def expectation(seq, nchances, gcprob):
    logexp = string_logprob(seq, gcprob) + numpy.log(nchances)
    return numpy.exp(logexp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        n = int(lines[0])
        seq = lines[1]
        gcs = [float(item) for item in lines[2].split()]
        nchances = n - len(seq) + 1

        expectations = [expectation(seq, nchances, gc) for gc in gcs]
        print(" ".join([str(e) for e in expectations]))
