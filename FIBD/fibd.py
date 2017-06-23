#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys

def mortal_fib(ngen, lifespan):
    pop_structure = [0] * lifespan
    pop_structure[0] = 1
    for i in range(1, ngen):
        offspring = sum(pop_structure[1:])
        pop_structure[1:] = pop_structure[0:-1]
        pop_structure[0] = offspring
    return sum(pop_structure)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        lines = f.readlines()
        line = "".join([l.strip() for l in lines])

        items = line.split()
        n, m = int(items[0]), int(items[1])

        result = mortal_fib(n, m)
        print(result)
