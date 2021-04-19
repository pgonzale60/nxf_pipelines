#!/usr/bin/env python3
"""
Filter pacbio reads with contigous repeat at their start or end and trim telomeric repeats.
This script relies on the high fidelity and strandedness of PacBio HiFi reads.
The strongest assumption is that the forward motif will only occur on the right
end of the reads, while the reverse complement of the motif can only
occur on the left side.

Usage:
        filter_telomeric_reads.py [--motif STR] [--times INT]
                                  [--out FILE] [--lacking FILE]

options:
    -m STR, --motif STR     telomeric repeat
                            [Default: TTAGGC]
    --times INT             minimum number of contiguous occurrences.
                            [Default: 50]
    -o FILE, --out FILE     filename for gzip compressed for telomeric reads.
                            [Default: telomericReads.fasta.gz]
    -l FILE, --lacking FILE filename for gzip compressed for non-telomeric reads.
"""

import gzip
import sys
from itertools import groupby
from docopt import docopt
from subprocess import Popen, PIPE
import shlex
import re

__author__ = "Pablo Manuel Gonzalez de la Rosa"
__version__ = '1.0.0'

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

def get_start_and_end_of_sequence(seq, size, nmer_size):
    return seq[0:size + nmer_size - 1] + ("N" * size) + seq[-(size + nmer_size - 1)::]

def get_sequence_start(seq, size):
    return seq[0:size]

def get_sequence_end(seq, size):
    return seq[-size::]

def self_concatenate_n_times(motif, n):
    return motif * n

def reverse_complement_sequence(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def defunct_trim_sequence_from_start(seq, motif):
    motif_size = len(motif)
    rev_motif = reversed(motif)
    motif_pos = 0
    updated_min_search_pos = 0 
    # expanding the search space up to almost twice the motif size
    # leaves space for inexact beggining of match,
    # which can correspond to degenerate sequence or 
    # incomplete motif at the beggining of sequence
    # allowing up to a 3 nucleotides insertion and
    # some motif mismatch
    updated_top_search_pos = motif_size + (motif_size - 1)
    while motif_pos >= 0:
        min_search_pos = updated_min_search_pos
        top_search_pos = updated_top_search_pos
        motif_pos = seq[min_search_pos:top_search_pos].find(motif)
        if motif_pos == 0:
            motif_pos = motif_size
        updated_min_search_pos += motif_pos
        updated_top_search_pos += motif_pos
    return seq[min_search_pos:]


def trim_sequence_from_start(seq, motif, min_occur):
    motif_size = len(motif)
    trimmedSeq = re.sub("(" + motif
                            + ".?)*.*("
                            + motif
                            + ".?){2,}"
                            + motif, "", seq)
    return trimmedSeq

def double_trim_sequence_from_start(seq, motif, min_occur):
    trimmed1 = trim_sequence_from_start(seq, motif, min_occur)
    doubleTrimmed = trim_sequence_from_start(trimmed1, motif, min_occur)
    return doubleTrimmed



if __name__ == "__main__":
    minLen              = 500
    args                = docopt(__doc__)
    outfile             = args['--out']
    nontelomfile        = args['--lacking']
    motif               = args['--motif']
    rev_motif           = reverse_complement_sequence(motif)
    motif_size          = len(motif)
    min_occur           = int(args['--times'])
    supermotif_size     = motif_size * int(args['--times'])
    supermotif          = motif      * int(args['--times'])
    rev_supermotif      = reverse_complement_sequence(supermotif)
    searchSpace         = minLen * 2
    write_non_telomeric = bool(nontelomfile)
    non_telomeric_reads = ''
    telomeric_reads     = ''

    for name, seq, qual in readfq(sys.stdin):
        read_size = len(seq)
        trimmed_sequence = ''
        if read_size >= searchSpace:
            seq_start = seq[0:searchSpace]
            seq_end   = seq[-searchSpace:]
            matches_start = re.finditer("("+rev_motif  + "){" + str(min_occur) +",}" , str(seq_start), re.I)
            matches_end = re.finditer("("+ motif  + "){" + str(min_occur) +",}" , str(seq_end), re.I)
            if any(matches_start):
                if not any(matches_end):
                    trimmed_sequence = trim_sequence_from_start(seq,
                                                                rev_motif, min_occur)
            elif any(matches_end):
                trimmed_sequence = trim_sequence_from_start(
                    reverse_complement_sequence(seq),
                 rev_motif, min_occur)
            if len(trimmed_sequence) > minLen:
                telomeric_reads += ">%s\n" % name
                telomeric_reads += "%s\n" % trimmed_sequence
            elif write_non_telomeric:
                non_telomeric_reads += ">%s\n" % name
                non_telomeric_reads += "%s\n" % seq

                
    with gzip.open(outfile, 'wt') as ofh:
        ofh.writelines(telomeric_reads)
    
    if write_non_telomeric:
        with gzip.open(nontelomfile, 'wt') as nofh:
            nofh.writelines(non_telomeric_reads)
