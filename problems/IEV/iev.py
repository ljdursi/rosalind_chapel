#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import numpy

def allele_count_probabilities(parent1, parent2):
    # calculate transition matrices
    T0 = numpy.zeros((3,3))
    T2 = numpy.zeros((3,3))

    for i in range(3):
        for j in range(3):
            T0[i, j] = (2 - i)*(2 - j)/4.0
            T2[i, j] = i*j/4.0

    allele_probs = numpy.zeros(3)
    allele_probs[0] = numpy.dot( parent1, numpy.dot(T0, parent2) )
    allele_probs[2] = numpy.dot( parent1, numpy.dot(T2, parent2) )
    allele_probs[1] = 1. - allele_probs[0] - allele_probs[2]
    return allele_probs

def frac_expected_dominant(parent1, parent2):
    probs = allele_count_probabilities(parent1, parent2)
    return probs[2] + probs[1]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-n', '--off', type=int, default=2)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        counts = [int(i) for i in lines[0].split()]

        AA = numpy.array([0., 0., 1.])
        Aa = numpy.array([0., 1., 0.])
        aa = numpy.array([1., 0., 0.])

        expectation = 0.0
        expectation += args.off * counts[0] * frac_expected_dominant(AA, AA)
        expectation += args.off * counts[1] * frac_expected_dominant(AA, Aa)
        expectation += args.off * counts[2] * frac_expected_dominant(AA, aa)
        expectation += args.off * counts[3] * frac_expected_dominant(Aa, Aa)
        expectation += args.off * counts[4] * frac_expected_dominant(Aa, aa)
        expectation += args.off * counts[5] * frac_expected_dominant(aa, aa)

        print(expectation)
