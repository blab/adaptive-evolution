from Bio import SeqIO
from Bio.Seq import Seq
from sys import stderr, stdin
import os
import argparse
#take in genome and wite out full sequence
#take in config file and alignment file as arguments
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--alignment", metavar="FASTA", default="results/aligned.fasta",
        help="The parse alignment output")
    parser.add_argument("--reference", metavar="GB", default="config/reference.gb",
        help="The reference genbank file")
    parser.add_argument('--percentage', metavar='N', type=float, nargs='+', default=0.8,
                    help='min percentage sequenced')
    parser.add_argument('--gene', metavar='N',nargs='+', default='all_genes',
                    help='gene sequenced')
    parser.add_argument('--output', metavar='N',nargs='+', default='result.fasta',
                    help='output file to store sequence')
    args = parser.parse_args()

output_files = {}
if args.gene[0].lower() != 'genome':
    # extract all keys and locations from genbank file
    genes = {}
    for seq_record in SeqIO.parse(args.reference, "genbank"):
        for feature in seq_record.features:
            if feature.type == 'CDS':
                if 'gene' in feature.qualifiers.keys():
                    gene_name = feature.qualifiers['gene']
                    gene_location = feature.location
                    genes[gene_name[0].lower()] = gene_location

    #if requested gene exists, only keep requested gene
    if args.gene[0].lower() in genes.keys():
        value = genes[args.gene[0].lower()]
        genes.clear()
        genes[args.gene[0].lower()] = value

    #generate output files in results/gene/aligned.fasta structure
    if len(genes.keys()) > 1:
        general_file_path = "results"
        for gene in genes.keys():
             dir_name = general_file_path + "/" + gene
             filename = dir_name + "/aligned.fasta"
             if not os.path.exists(dir_name):
                 os.makedirs(dir_name)
             output_files[gene] = open(filename, "w+")
    else:
        if not os.path.exists(os.path.dirname(args.output[0])):
            os.makedirs(os.path.dirname(args.output[0]))
        output_files[list(genes.keys())[0]] = open(args.output[0], "w+")

    #sequence dictionary to store sequence strings
    sequences = {}

    #parsing and concatenating sequences, does not write strands that are sequenced below given percentage
    fasta_sequences = SeqIO.parse(open(args.alignment),'fasta')
    for record in fasta_sequences:
        for gene in genes.keys():
            if not gene in sequences:
                sequences[gene]  = ""
            seq = genes[gene].extract(record.seq)
            #calculate if percent of Ns is less than percentage
            if (seq.lower().count('n')/len([x for x in seq if x.isalpha()]) <= args.percentage[0]):
                new_gene = sequences[gene]  + ">" + record.id + "\n" + seq + "\n"
                sequences[gene] = new_gene

    #write sequences to output files
    for file in output_files.keys():
        output_files[file].write(str(sequences[file]))
else:
    #make file for full length genome
    if not os.path.exists(os.path.dirname(args.output[0])):
        os.makedirs(os.path.dirname(args.output[0]))
    output_file = open(args.output[0], "w+")
    #parsing and concatenating sequences, does not write strands that are sequenced below given percentage
    long_seq = ""
    seq = ""
    fasta_sequences = SeqIO.parse(open(args.alignment),'fasta')
    for record in fasta_sequences:
        seq = record.seq
        #print(seq)
        #calculate if percent of Ns is less than percentage
        if (seq.lower().count('n')/len([x for x in seq if x.isalpha()]) <= args.percentage[0]):
            long_seq = long_seq  + ">" + record.id + "\n" + seq + "\n"
    #write to output file
    output_file.write(str(long_seq))
