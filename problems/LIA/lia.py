#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import numpy
from scipy.stats import binom


def allele_count_probabilities(ngen, first_parent, other_parents):
    # calculate transition matrices
    T0 = numpy.zeros((3, 3))
    T2 = numpy.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            T0[i, j] = (2 - i)*(2 - j)/4.0
            T2[i, j] = i*j/4.0

    # for first_parent = other_parents = [0,1,0],
    # result will always be [1/4, 1/2, 1/4]
    allele_probs = first_parent[:]
    for gen in range(ngen):
        curprobs = allele_probs[:]
        allele_probs[0] = numpy.dot(other_parents, numpy.dot(T0, curprobs))
        allele_probs[2] = numpy.dot(other_parents, numpy.dot(T2, curprobs))
        allele_probs[1] = 1. - allele_probs[0] - allele_probs[2]
    return allele_probs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-n', '--nalleles', type=int, default=2)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        line = "".join([l.strip() for l in lines])

        items = line.split()
        k, n = int(items[0]), int(items[1])

        tom = numpy.array([0., 1., 0.])
        others = numpy.array([0., 1., 0.])
        probabilities = allele_count_probabilities(k, tom, others)
        prob_one_allele = probabilities[1]**args.nalleles
        n_exactly_one = 1 - binom.cdf(n-1, 2**k, prob_one_allele)
        print(n_exactly_one)
