#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys


def multiset_from_string(line):
    items = line.split()
    nums = [int(item) for item in items if not item == ""]
    nums.sort()
    return nums


def infer_set_from_diff_multiset(dmultiset):
    values = [0, dmultiset[-1]]
    diffs = dmultiset[:-1]

    def update_diffs_with_new_value(alldiffs, vals, newval):
        newdiffs = [abs(v - newval) for v in vals]
        alldiffs_copy = alldiffs[:]
        for diff in newdiffs:
            try:
                alldiffs_copy.remove(diff)
            except:
                return None
        return alldiffs_copy

    def solve(values, diffs):
        if diffs is None:
            return None

        if len(diffs) == 0:
            return values

        nextdiff = diffs[-1]
        leftval = values[-1] - nextdiff
        rightval = values[0] + nextdiff

        leftdiffs = update_diffs_with_new_value(diffs, values, leftval)
        leftvals = solve(sorted(values + [leftval]), leftdiffs)
        if leftvals is not None:
            return leftvals

        rightdiffs = update_diffs_with_new_value(diffs, values, rightval)
        rightvals = solve(sorted(values + [rightval]), rightdiffs)
        if rightvals is not None:
            return rightvals

        return None

    return solve(values, diffs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for line in f.readlines():
            diffs = multiset_from_string(line.strip())
            vals = infer_set_from_diff_multiset(diffs)
            print(" ".join([str(v) for v in vals]))
