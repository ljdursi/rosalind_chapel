#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys


def complement(base):
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return complements[base]


def reversecomplement(line):
    rc = ''.join([complement(c) for c in reversed(line)])
    return rc


def fasta_iterator(infile):
    started = False
    label = ""
    sequence = ""

    for line in infile:
        line = line.strip()
        if line[0] == '>':
            if started:
                yield label, sequence
            started = True
            label = line[1:]
            sequence = ""
        else:
            sequence += line

    if started:
        yield label, sequence


def growpalindrome(seq, pivotoffset, knownlen, maxlen):
    """
    Given a known complement palindrome seed at a given
    location in a sequence, attempt to find the full extent
    of the palindrome up to length maxlen.
    """
    palindromelen = knownlen
    halflen = knownlen // 2

    # trust, but verify
    firsthalf = seq[pivotoffset-halflen+1:pivotoffset+1]
    secondhalf = seq[pivotoffset+1:pivotoffset+1+halflen]
    assert firsthalf == reversecomplement(secondhalf)

    for l in range(halflen+1, maxlen//2+1):
        first, last = pivotoffset-l+1, pivotoffset+l
        if first < 0 or last >= len(seq):
            break
        if not seq[first] == complement(seq[last]):
            break
        palindromelen += 2

    return palindromelen


def findallminpalindromes(seq, minlen):
    halflen = minlen//2

    minpalindromes = []
    for i in range(len(seq)-minlen+1):
        kmer = seq[i:i+halflen]
        nextkmer = seq[i+halflen:i+minlen]
        if reversecomplement(nextkmer) == kmer:
            minpalindromes.append(i)

    return minpalindromes


def allpalindromes(seq, minlen, maxlen):
    palindromes = []
    minpalindromes = findallminpalindromes(seq, minlen)
    for minpalindrome in minpalindromes:
        pivotoffset = minpalindrome + minlen//2 - 1
        pmaxlen = growpalindrome(seq, pivotoffset, minlen, maxlen)
        for plen in range(minlen, pmaxlen+1, 2):
            newstart = minpalindrome-(plen-minlen)//2
            palindromes.append((newstart, plen, seq[newstart:newstart+plen]))
    return palindromes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--minlen', type=int, default=4)
    parser.add_argument('-M', '--maxlen', type=int, default=12)
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    args = parser.parse_args()

    with args.infile as f:
        for label, sequence in fasta_iterator(f):
            res = allpalindromes(sequence, args.minlen, args.maxlen)
            res.sort(key=lambda x: x[0])
            for item in res:
                print(item[0]+1, item[1])
