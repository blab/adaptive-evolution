from augur.utils import read_node_data
import pandas as pd
from Bio import SeqIO
import re

#This just adds the correct sequence length, without - or N, for those that were actually
#used in the analysis (passed 'filter') to the metadata

if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='add seq-len',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--seqs-in', help="input aligned sequences")           
    parser.add_argument('--meta-in', help="input meta file")
    parser.add_argument('--meta-out', help="output metadata with seq-lens based on alignment and with N and - removed")
    args = parser.parse_args()

    orig_meta = args.meta_in #"results/metadata-ages.tsv"
    seqs = args.seqs_in #"results/aligned_vp1.fasta"
    meta = pd.read_csv(orig_meta, sep='\t', index_col=False)

    record_dict = SeqIO.to_dict(SeqIO.parse(seqs, "fasta"))

    #only re-do sequence length for sequences that passed 'filter' (not all in meta may have)
    for i, row in meta.iterrows():
        if row.strain in record_dict:
            str_seq = str(record_dict[row.strain].seq)
            meta.loc[meta.strain == row.strain, 'seq-len'] = len(re.sub("-|N","", str_seq))

    meta.to_csv(args.meta_out, sep="\t", index=False)
