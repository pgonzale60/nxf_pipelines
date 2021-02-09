#!/usr/bin/env python3
"""Update coordinates and remove seq name suffix from chunked blast results."""

import re
from collections import defaultdict
from docopt import docopt

docs = """
Unchunk BLAST results.

Usage: ./unchunk_blast.py [--count INT] [--in TSV] [--out TSV]

Options:
    --in TSV     input filename.
    --out TSV    output filename.
    --count INT  number of results to keep per chunk. [Default: 10]
"""

__author__ = "Richard Challis, Pablo Manuel Gonzalez de la Rosa"



if __name__ == '__main__':
    args = docopt(docs)
    try:
        lines = defaultdict(dict)
        chunk_counts = defaultdict(int)
        with open(args['--in'], 'r') as fh:
            for line in fh.readlines():
                if '\n' not in line:
                    line += '\n'
                fields = line.split('\t')
                if fields[0]:
                    name, start = re.split('_-_', fields[0])
                    fields[0] = name
                    fields[1] = str(int(fields[1])+int(start))
                    fields[2] = str(int(fields[2])+int(start))
                    if start not in lines[name]:
                        lines[name][start] = []
                        chunk_counts[name] += 1
                    lines[name][start].append('\t'.join(fields))
        with open(args['--out'], 'w') as ofh:
            for name in chunk_counts.keys():
                length = len(lines[name])
                n = int(args['--count'])
                for start in lines[name].keys():
                    for i in range(n):
                        if i < len(lines[name][start]):
                            ofh.write("%s" % (lines[name][start][i]))
    except Exception as err:
        logger.error(err)
        exit(1)
