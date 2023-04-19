
from augur.utils import read_node_data
import pandas as pd
from Bio import SeqIO

#after this may want to save a translated (AA) copy. lock all gaps (-)
#can get into sequential format in R:
#seqs <- read.fasta("C:/Users/Emma/bash/github/enterovirus_vp1_test/results/seqs-subgeno-trans.fasta", forceDNAtolower=F, as.string=T)
#write.fasta(seqs, attr(seqs,"name"), "C:/Users/Emma/bash/github/enterovirus_vp1_test/results/seqs-subgeno-trans2.fasta", nbchar=6000)


if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='add subgenogroup',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--seqs-in', help="input sequences")           
    parser.add_argument('--meta-in', help="input meta file")
    parser.add_argument('--clades', help="clades JSON file")
    parser.add_argument('--meta-out', help="output metadata just with subgenogroup added")
    parser.add_argument('--seqs-out', help="output sequences just with subgenogroup added")
    args = parser.parse_args()

    orig_meta = args.meta_in #"results/metadata-ages.tsv"
    clade_info = args.clades #"results/clades_vp1.json"
    seqs = args.seqs_in #"results/aligned_vp1.fasta"

    meta = pd.read_csv(orig_meta, sep='\t', index_col=False)

    clade_node = read_node_data(clade_info)
    clade_node = clade_node["nodes"]

    record_dict = SeqIO.to_dict(SeqIO.parse(seqs, "fasta"))

    to_exclude = []

    for i, row in meta.iterrows():
        #only do this if the strain is in the clades
        if row.strain in clade_node.keys() and row.strain in record_dict.keys():
            meta.loc[meta.strain == row.strain, 'subgenogroup'] = clade_node[row.strain]['clade_membership']
            if clade_node[row.strain]['clade_membership'] == 'unassigned':
                meta.loc[meta.strain == row.strain, 'subgenogroup'] = ""
            else:
                record_dict[row.strain].id = record_dict[row.strain].name + "-Clade" + clade_node[row.strain]['clade_membership']
                record_dict[row.strain].name = "" #record_dict[row.strain].id
                record_dict[row.strain].description = ""
        else:
            to_exclude.append(row.strain)
            if row.strain in record_dict.keys():
                del record_dict[row.strain]


    seq_values = [v for v in record_dict.values()]
    #SeqIO.write(seq_values, "results/seqs-subgeno.fasta", "fasta")
    SeqIO.write(seq_values, args.seqs_out, "fasta")

    mask = meta['strain'].isin(to_exclude)
    #meta[~mask].to_csv("results/metadata-subgeno.tsv", sep='\t', index=False)
    meta[~mask].to_csv(args.meta_out, sep='\t', index=False)
