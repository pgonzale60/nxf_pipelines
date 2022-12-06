#!/usr/bin/env python3
"""
Identify softmasked positions from fasta and write in BED format.

Usage:
        soft_mask2bed.py [--fasta FILE] [--out FILE]

options:
    -f FILE, --fasta FILE   filename of softmasked fasta file
    -o FILE, --out FILE     filename for softmasked regions in BED format.
                            [Default: softmasked.bed]
"""

import sys
from docopt import docopt

# This file will generate a bedfile of the masked regions a fasta file.

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

#via https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
def list2groups(L):
    if len(L) < 1:
        return
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last: # Part of the group, bump the end
            last = n
        else: # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last # Yield the last group

def masked2bed(infasta, outbed):
    masked = []
    maskedSize = 0
    with open(infasta, 'rU') as infile:
        for header, seq, qual in readfq(infile):
            # via https://github.com/nextgenusfs/redmask/blob/master/redmask.py
            if ' ' in header:
                ID = header.split(' ')[0]
            else:
                ID = header
            for i,c in enumerate(seq):
                if c.islower():
                    masked.append(i) #0 based
                    maskedSize += 1
                
    if maskedSize > 0: # only if softmasked
        with open(outbed, 'w') as outfile:
            repeats = list(list2groups(masked))
            for item in repeats:
                if len(item) == 2:
                    outfile.write('{:}\t{:}\t{:}\n'.format(ID, item[0], item[1]))


if __name__ == "__main__":
    args      = docopt(__doc__)
    fastafile = args['--fasta']
    outfile   = args['--out']

    masked2bed(fastafile, outfile)
