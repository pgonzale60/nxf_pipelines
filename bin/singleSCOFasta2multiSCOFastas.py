#!/usr/bin/env python3
import os
import shutil
import argparse
import pyfastx
import itertools
from collections import Counter
import pandas as pd


def parse_fasta(fasta_file):
    busco_seqs = pyfastx.Fasta(fasta_file)
    ids = busco_seqs.keys()
    busco_names = pd.DataFrame(data=list(ids), columns=['seqNames'])
    busco_names = pd.DataFrame(busco_names.seqNames.str.split("_", 1).tolist(),
                    columns=['buscoId', 'sampleId'])
    return busco_names, busco_seqs

def count_single_copy_buscos(busco_seqs):
    scb_counts = (
            busco_seqs['buscoId']
            .groupby(busco_seqs['buscoId'])
            .transform('size')
        )
    return scb_counts

def get_usable_buscos(busco_seqs, proportion):
    n_taxa = len(busco_seqs.sampleId.unique())
    min_taxon_count = int(round(n_taxa * proportion, 0))
    print("[STATUS] Selecting BUSCOs present and single-copy in ≥ " + str(min_taxon_count) + " taxa...")
    scb_count = count_single_copy_buscos(busco_seqs)
    usable_scb = busco_seqs[scb_count > min_taxon_count]
    usable_scb.loc[:, ('seqNames')] = usable_scb.apply(lambda x: '_'.join(x), axis=1)
    return usable_scb


def create_output_fastas(outdir, usable_scb, sco_seqs, suffix):
    if os.path.exists(outdir):
        print("[WARING] Removing existing output directory ('" + outdir + "') and its contents...")
        shutil.rmtree(outdir)
    os.mkdir(outdir)
    selectedBuscos = usable_scb.groupby('buscoId')
    total, count = selectedBuscos.ngroups, 0
    print("[STATUS] Writing " + str(total) + " " + suffix + " FASTA files to output directory ('" + (outdir) + "'). This may take a moment...")
    for group_name, df_group in selectedBuscos:
        with open(outdir + "/" + group_name + suffix, 'w') as outfile:
            for row_index, row in df_group.iterrows():
                outfile.write(">" + row['sampleId'] + "\n")
                outfile.write(str(sco_seqs[row['seqNames']]) + "\n")
        print("[STATUS] \t" + str(count) + "/" + str(total) + " files written.", end='\r')
        count += 1
    print("[STATUS] \t" + str(count) + "/" + str(total) + " files written.")



if __name__ == "__main__":
    SCRIPT = "singleSCOFasta2multiSCOFastas.py"
    ### argument set up 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ortho_seqs", type=str, help = "directory containing a set BUSCO results directories (e.g. one per taxon)", required=True)
    parser.add_argument("-o", "--outdir", type=str, help = "output directory for FASTAs (default: b2f_output)", default="b2f_output")
    parser.add_argument("-s", "--seqtype", type=str, help = "your chosen sequence type (default: protein)", choices={"protein", "nucleotide"}, default="protein")
    parser.add_argument("-p", "--proportion", type=float, help = "proportion of taxa required for a given BUSCO to be output as FASTA (default: 1.0)", default=1.0)
    args = parser.parse_args()
    seqs_file = args.ortho_seqs
    seqtype = args.seqtype
    proportion = args.proportion
    outdir = args.outdir
    # get suffix based on seqtype
    if seqtype == 'protein':
        suffix = ".faa"
    else:
        suffix = ".fna"
    # proportion
    print("singleSCOFasta2multiSCOFastas.py parameters:")
    print("\tinput_file: " + seqs_file)
    print("\toutput_dir: " + outdir)
    print("\tseqtype: " + seqtype)
    print("\tproportion: " + str(proportion))
    print("")
    ### run functions
    sco_names, sco_seqs = parse_fasta(seqs_file)
    usable_scb = get_usable_buscos(sco_names, proportion)
    create_output_fastas(outdir, usable_scb, sco_seqs, suffix)
    print("[STATUS] singleSCOFasta2multiSCOFastas.py completed successfully.")
