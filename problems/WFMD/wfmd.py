#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import numpy
from scipy.special import binom


def transition_matrix(npop):
    twoN = 2*npop
    T = numpy.zeros((twoN+1, twoN+1))

    twon_choose_j = numpy.zeros(twoN+1)
    for j in range(twoN+1):
        twon_choose_j[j] = binom(twoN, j)

    for i in range(twoN+1):
        for j in range(twoN+1):
            T[j, i] = twon_choose_j[j]*(1.0*i/twoN)**j
            T[j, i] *= (1.0*(twoN-i)/twoN)**(twoN-j)

    return T


def evolve(initial_vec, npop, ngen):
    v = initial_vec[:]
    T = transition_matrix(npop)
    for _ in range(ngen):
        v = numpy.dot(T, v)

    return v


def initial_vector(npop, m):
    v = numpy.zeros((2*npop+1))
    v[m] = 1.0
    return v


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            npop, m, g, k = (int(item) for item in line.split())
            v = initial_vector(npop, m)
            v = evolve(v, npop, g)
            print(numpy.sum(v[:-k]))
