#! /usr/local/Caskroom/mambaforge/base/envs/nextstrain/bin/python

# Get all insertions and code as binary vector: Use metadata for this
# Load tree
# Reconstruct using tree time

import json
from collections import defaultdict

import augur.ancestral as ancestral
import Bio.Align
import click
import pandas as pd
from BCBio import GFF
from mergedeep import merge


@click.command()
@click.option("--insertions", type=str)
@click.option("--tree", type=str)
@click.option("--output", type=click.File("w"))
@click.option("--genemap", default = "config/genemap.gff", type=click.File("r"))
def main(insertions, tree, output, genemap):
    meta = pd.read_csv(insertions, sep=",", index_col=0)
    # Find all unique insertions, create list

    result = {}

    def split_csv_field(x):
        """Split csv field into list, with empty strings yielding []"""
        return x.split(";") if type(x) == str else []

    # Construct of insertions
    characters = {
        character for characters in meta['insertions'] for character in split_csv_field(characters)
    }

    # Dictionary that maps property to position in pseudosequence
    mapping = {}
    for pos, character in enumerate(characters):
        mapping[character] = pos

    inverse_mapping = {v: k for k, v in mapping.items()}

    # Take list of characters and turn into binary vector
    # Convention: A = 0, G = 1
    def character_list_to_vector(character_list_csv):
        """Turn comma separated list into binary vector"""
        character_list = split_csv_field(character_list_csv)
        vector = ["A"] * len(characters)
        for character in character_list:
            vector[mapping[character]] = "G"
        return "".join(vector)

    meta = meta.assign(trait_vector=meta['insertions'].apply(character_list_to_vector))


    # Print mapping
    # print(mapping)
    # Print all sequences with insertions
    # print(meta.insertions_vector[meta.insertions_vector.apply(lambda x: "G" in x) > 0])

    # Transform binary vector into BioSeq alignment
    alignment = Bio.Align.MultipleSeqAlignment(
        [
            Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(vector), id=strain, name=strain)
            for strain, vector in meta.trait_vector.iteritems()
        ]
    )
    # print(alignment)

    from treetime import TreeAnc

    tt = TreeAnc(tree=tree, aln=alignment)
    tt.infer_ancestral_sequences()
    nodes_with_mutations = ancestral.collect_mutations_and_sequences(tt)

    gm = list(GFF.parse(genemap))[0]
    def position_to_gene(position):
        """Convert position to gene using genemap"""
        gene_name = codon_number = reading_frame = None
        for feature in gm.features:
            if position in feature.location:
                gene_name = feature.qualifiers["gene_name"][0]
                codon_number = (position - feature.location.start) // 3 + 1
                reading_frame = (position - feature.location.start - 1) % 3
        return {"gene_name": gene_name, "codon_number": codon_number, "reading_frame": reading_frame}

    def mut_to_str(mut: str) -> dict:
        """Convert mutation to string"""

        indel_type = "ins" if mut[0] == "A" else "rev ins"
        insertion = inverse_mapping[int(mut[1:-1]) - 1].split(":")
        insertion_position = insertion[0]
        inserted_nucleotides = insertion[1]
        gene_position = position_to_gene(int(insertion_position))
        if gene_position["gene_name"] is None:
            position_string = ""
        else:
            frame_string = f"/{gene_position['reading_frame']}" if gene_position["reading_frame"] != 0 else ""
            position_string = f" ({gene_position['gene_name']}:{gene_position['codon_number']}{frame_string})"
        return (indel_type, insertion_position + inserted_nucleotides + position_string)


    for node in nodes_with_mutations["nodes"].values():
        list_of_pairs = list(map(mut_to_str, node["muts"]))
        node["muts"] = defaultdict(list)
        for a,b in list_of_pairs:
            node["muts"][a].append(b)
        node["muts"] = dict(node["muts"])

    # for key, value in nodes_with_mutations["nodes"].items():
    #     if value["muts"] != {}:
    #         print(f"{key}: {value['muts']}")

    merge(result, nodes_with_mutations)

    json.dump(result , output)


if __name__ == "__main__":
    main()
