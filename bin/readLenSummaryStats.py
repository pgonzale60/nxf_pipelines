#!/usr/bin/env python3

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import groupby 
import matplotlib
import argparse
import gzip

matplotlib.use('AGG') # Save matplots as PNGs


def parse_fasta(fasta_file):
    stpmx = []
    with gzip.open(fasta_file) as handle:
        for header, group in groupby(handle, lambda x:x.startswith(b'>')):
            if not header:
                stpmx.append(sum(map(lambda x: len(x.strip()), group)))
    return np.array(stpmx)

def length_stats(ndLengths):
    seqYield = np.sum(ndLengths)
    nSeqs = ndLengths.size
    [n50, n95, n100] = np.quantile(ndLengths, [.5,.95, 1])
    return seqYield, nSeqs, n50, n95, n100

def plot_length_dist(ndLengths, outfile, n50, n95, n100):
    plt.ioff()
    ax = sns.histplot(data = ndLengths, fill=False, element="poly")
    kde_x, kde_y = ax.lines[0].get_data()
    p1 = plt.axvline(x=n50,color='#62a8dd')
    p2 = plt.axvline(x=n95,color='#0d5698')
    dot1 = ax.fill_between(kde_x, kde_y, where=(kde_x<=n100), 
                    interpolate=True, color='#0d5698')
    dot2 = ax.fill_between(kde_x, kde_y, where=(kde_x<n95) , 
                    interpolate=True, color='#62a8dd')
    dot3 = ax.fill_between(kde_x, kde_y, where=(kde_x<n50), 
                    interpolate=True, color='#9fd6f6')
    ax.legend([dot2, dot3], ['>= N50', '>= 95th quantile'])
    ax.set(xlabel='Read length (bp)')
    ax.figure.savefig(outfile)

if __name__ == "__main__":
    SCRIPT = "readLenSummaryStats.py"
    ### argument set up 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fasta", type=str, help = "fasta file to summarize", required=True)
    parser.add_argument("-o", "--outfile", type=str, help = "name for plot file", required=True)
    args = parser.parse_args()
    seqs_file = args.fasta
    outfile = args.outfile
    ### run functions
    readLen = parse_fasta(seqs_file)
    seqYield, nSeqs, n50, n95, n100 = length_stats(readLen)
    plot_length_dist(readLen, outfile, n50, n95, n100 )
    print(*[outfile, seqYield, nSeqs, n50, n95, n100], sep='\t')
