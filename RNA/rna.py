#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
from collections import defaultdict

def transcribe(line):
    result = line.replace('T', 'U')
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        line = "".join([l.strip() for l in lines])

        result = transcribe(line)
        print(result)
