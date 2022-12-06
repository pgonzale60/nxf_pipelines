#!/usr/bin/env python3
"""
Split reads into files containing X number of reads.

Usage:
        cat reads | split_reads.py [-s INT]

options:
    -s INT, --size INT     Number of reads
                            [Default: 50000]
"""

import sys
from docopt import docopt

__author__ = "Pablo Manuel Gonzalez de la Rosa"
__version__ = '1.1.1'

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





if __name__ == "__main__":
    args    = docopt(__doc__)
    outSize = args['--size']
    outseq  = ''

    for name, seq, qual in readfq(sys.stdin):
        read_size = len(seq)
        trimmed_sequence = ''
        if read_size >= searchSpace:
            matches_start = seq.find(rev_motif, 0, motif_size * 3)
            matches_end = seq.find(motif, -motif_size * 3)
            if matches_start >= 0:
                if not matches_end >= 0:
                    telomere_pos = seq.rfind(rev_motif * min_occur) + (motif_size * min_occur)
                    trimmed_sequence = seq[telomere_pos:]
            elif matches_end >= 0:
                telomere_pos = seq.find(motif * min_occur)
                trimmed_sequence = reverse_complement_sequence(seq[0:telomere_pos])
            if len(trimmed_sequence) > minLen:
                telomeric_reads += ">%s\n" % name
                telomeric_reads += "%s\n" % trimmed_sequence
            elif write_non_telomeric:
                non_telomeric_reads += ">%s\n" % name
                non_telomeric_reads += "%s\n" % seq

                
    with open(outfile, 'w') as ofh:
        ofh.writelines(telomeric_reads)
    
    if write_non_telomeric:
        with open(nontelomfile, 'w') as nofh:
            nofh.writelines(non_telomeric_reads)
