#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys


def differences(a, b, delta=1.e-6):
    n = len(a)
    m = len(b)
    diffs = [0.0]*(n*m)

    i = 0
    for aa in a:
        for bb in b:
            diffs[i] = aa-bb
            i += 1

    diffs.sort()

    last = -999999999999
    count = 0
    multiplicities = []
    for diff in diffs:
        if abs(diff-last) < delta:
            count += 1
        else:
            multiplicities.append((count, last))
            last = diff
            count = 1

    return multiplicities[1:] + [(count, last)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-d', '--delta', default=1.e-6)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines]

        a = [float(i) for i in lines[0].split()]
        b = [float(j) for j in lines[1].split()]

        diffs = differences(a, b, args.delta)
        diffs.sort()
        modediff = diffs[-1]
        print(modediff[0])
        print(abs(modediff[1]))
