#!/usr/bin/env python3
"""
Filter pacbio reads with contigous repeat at their start or end and trim telomeric repeats.
This script realies on the high fidelity and strandedness of PacBio HiFi reads.
The strongest assumption is that the forward motif will only occur on the right
end of the reads, while the reverse complement of the motif can only
occur on the left side.

Usage:
        filter_telomeric_reads.py [--in FASTA] [--motif STR] [--times INT]
                                  [--out FILE] [--lacking FILE]

options:
    -i FASTA, --in FASTA    input gzip compressed FASTA file.
    -m STR, --motif STR     telomeric repeat
                            [Default: TTAGGC]
    --times INT             minimum number of contiguous occurrences.
                            [Default: 50]
    -o FILE, --out FILE     filename for gzip compressed for telomeric reads.
                            [Default: telomericReads.fasta.gz]
    -l FILE, --lacking FILE filename for gzip compressed for non-telomeric reads.
                            [Default: nontelomericReads.fasta.gz]
"""

import gzip
from itertools import groupby
from docopt import docopt
from subprocess import Popen, PIPE
import shlex

__author__ = "Richard Challis, Pablo Manuel Gonzalez de la Rosa"
__version__ = '0.3.1'

def read_fasta(fastafile):
    """Read fasta"""
    cmd = "pigz -dc %s" % fastafile
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

def trim_sequence_from_start(seq, motif):
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



if __name__ == "__main__":
    args                = docopt(__doc__)
    outfile             = args['--out']
    nontelomfile        = args['--lacking']
    motif               = args['--motif']
    rev_motif           = reverse_complement_sequence(motif)
    motif_size          = len(motif)
    supermotif_size     = motif_size * int(args['--times'])
    supermotif          = motif      * int(args['--times'])
    rev_supermotif      = reverse_complement_sequence(supermotif)
    write_non_telomeric = bool(nontelomfile)
    non_telomeric_reads = ''
    telomeric_reads     = ''

    for seq in read_fasta(args['--in']):
        if len(seq['seq']) >= 2 * supermotif_size:
            seq_start = get_sequence_start(seq['seq'], supermotif_size * 2)
            seq_end = get_sequence_end(seq['seq'], supermotif_size * 2)
            begins_with_telomere = rev_supermotif in seq_start or supermotif in seq_start
            ends_with_telomere   = supermotif in seq_end or rev_supermotif in seq_end
            if begins_with_telomere:
                if not ends_with_telomere:
                    telomeric_reads += ">%s\n" % seq['title']
                    trimmed_sequence = trim_sequence_from_start(seq['seq'], rev_motif)
                    telomeric_reads += "%s\n" % trimmed_sequence
            elif ends_with_telomere:
                telomeric_reads += ">%s\n" % seq['title']
                trimmed_sequence = trim_sequence_from_start(
                    reverse_complement_sequence(seq['seq']),
                 rev_motif)
                telomeric_reads += "%s\n" % trimmed_sequence
            elif write_non_telomeric:
                non_telomeric_reads += ">%s\n" % seq['title']
                non_telomeric_reads += "%s\n" % seq['seq']

                
    with gzip.open(outfile, 'wt') as ofh:
        ofh.writelines(telomeric_reads)
    
    if write_non_telomeric:
        with gzip.open(nontelomfile, 'wt') as nofh:
            nofh.writelines(non_telomeric_reads)
