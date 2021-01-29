#!/usr/bin/env python3
"""
Filter sequences with contigous repeat at their start or end.

Usage:
        filter_telomeric_reads.py [--in FASTA] [--motif STR] [--times INT]
                                   [--out FILE]

options:
    -i FASTA, --in FASTA    input gzip compressed FASTA file.
    -m STR, --motif STR     telomeric repeat
                            [Default: TTAGGC]
    --times INT             minimum number of contiguous occurrences.
                            [Default: 50]
    -o FILE, --out FILE     output filename for gzip compressed FASTA.
                            [Default: telomericReads.fasta.gz]
"""

import pyfastx
import gzip
from itertools import groupby
from docopt import docopt
from subprocess import Popen, PIPE
import shlex

__author__ = "Richard Challis, Pablo Manuel Gonzalez de la Rosa"
__version__ = '0.1.0'

def read_fasta(fastafile):
    """Read fasta"""
    cmd = "zcat %s" % fastafile
    # TODO: read gzipped files if needed
    # cmd = "pigz -dc %s" % fastafile
    title = ''
    seq = ''
    with Popen(shlex.split(cmd), encoding='utf-8', stdout=PIPE, bufsize=4096) as proc:
        faiter = (x[1] for x in groupby(proc.stdout, lambda line: line[0] == '>'))
        for header in faiter:
            title = header.__next__()[1:].strip().split()[0]
            seq = ''.join(map(lambda s: s.strip(), faiter.__next__()))
    yield {'title': title, 'seq': seq}

def get_start_and_end_of_sequence(seq, size):
    return seq[0:size] + ("N" * size) + seq[-size::]

def self_concatenate_n_times(motif, n):
    return motif * n

def reverse_sequence(seq):
    return seq[::-1]





if __name__ == "__main__":
    args = docopt(__doc__)

    supermotif = self_concatenate_n_times(args['--motif'], int(args['--times']))
    rev_supermotif = reverse_sequence(supermotif)
    supermotif_size = len(supermotif)

    telomeric_reads = ''
    for seq in read_fasta(args['--in']):
        if len(seq['seq']) >= 2 * supermotif_size:
            start_and_end_seq = get_start_and_end_of_sequence(seq['seq'])
            if supermotif in start_and_end_seq or rev_supermotif in start_and_end_seq:
                telomeric_reads += ">%s\n" % seq['title']
                telomeric_reads += "%s\n" % seq['seq']
    
    outfile = args['--out']
    with gzip.open(outfile, 'w') as ofh:
        ofh.writelines(telomeric_reads)