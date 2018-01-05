#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import re


def set_from_string(line):
    items = re.split(',|{|}| +', line)
    nums = [int(item) for item in items if not item == ""]
    return set(nums)


def printset(inputset):
    output = "{"
    output += ", ".join([str(i) for i in inputset])
    output += "}"
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        lines = [l.strip() for l in lines]

        n = int(lines[0])
        a = set_from_string(lines[1])
        b = set_from_string(lines[2])

        u = set(range(1, n+1))
    
        printset(a.union(b))
        printset(a & b)
        printset(a - b)
        printset(b - a)
        printset(u - a)
        printset(u - b)
        
