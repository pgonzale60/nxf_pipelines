#!/usr/bin/env python3

"""
Filter pacbio reads with contigous repeat at their start or end and trim telomeric repeats.
This script relies on the high fidelity and strandedness of PacBio HiFi reads.
The strongest assumption is that the forward motif will only occur on the right
end of the reads, while the reverse complement of the motif can only
occur on the left side.

Usage:
        filter_telomeric_reads.py [--string STR] [--times INT]
                                  [--out FILE] [--buffered INT]
                                  [--min_len INT]

options:
    -s STR, --string STR    telomeric repeat
                            [Default: TTAGGC]
    --times INT             minimum number of contiguous occurrences.
                            [Default: 3]
    -m INT, --min_len INT   Minimum length of trimmed read
                            [Default: 100]
    -b INT, --buffered INT  Length of sequence to be kept in memory
                            [Default: 10000000]
    -o FILE, --out FILE     filename for telomeric reads.
                            [Default: telomericReads.fasta]
"""


from trim import readfq, reverse_complement_sequence, trim_seq

import sys
import os
from docopt import docopt

if __name__ == "__main__":
    args                = docopt(__doc__)
    motif               = args['--string']
    min_occur           = int(args['--times'])
    min_len             = int(args['--min_len'])
    buffer_size         = int(args['--buffered'])
    outfile             = args['--out']
    rev_motif           = reverse_complement_sequence(motif)
    longer_motif        = motif * min_occur
    longer_rev_motif    = rev_motif * min_occur
    motif_size          = len(motif)
    telomeric_reads     = ''
    seq_in_mem_len      = 0

    if os.path.exists(outfile):
        print("[WARNING] Overwriting existing file: {}".format(outfile))
        os.remove(outfile)

    for name, seq, qual in readfq(sys.stdin):
        trimmed_sequence = trim_seq(seq, min_occur, motif_size,
        motif, rev_motif, longer_motif, longer_rev_motif)
        trimmed_len = len(trimmed_sequence)
        if trimmed_len > min_len:
            telomeric_reads += ">%s\n" % name
            telomeric_reads += "%s\n" % trimmed_sequence
            seq_in_mem_len += trimmed_len
            if seq_in_mem_len > buffer_size:
                with open(outfile, 'a') as ofh:
                    ofh.writelines(telomeric_reads)
                telomeric_reads = ''
                seq_in_mem_len = 0
    
    if seq_in_mem_len > 0:
        with open(outfile, 'a') as ofh:
                    ofh.writelines(telomeric_reads)