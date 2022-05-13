import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from datetime import timedelta, datetime


def get_parent(tree, child_clade):
    """
    Function that returns the path from root to specified clade
    """
    node_path = tree.get_path(child_clade)
    return node_path


def nuc_changes_from_root(muts_on_path):
    """
    From all the of the nucleotide changes that have occurred on the path from root to branch,
    find the most recent nuc mutation at each site (giving the genotype at the branch)
    """

    final_muts_from_root = {}

    # overwrites genotypes at pos in historical order
    for x in muts_on_path:
        x_pos = int(x[1:-1])
        final_muts_from_root[x_pos] = x[-1]


    return final_muts_from_root


def determine_synonymous(nuc_muts_on_branch, parent_diffs_from_ref, reference_gene_locations, reference_gene_codon, reference_sequence_nt, reference_sequence_aa):
    """
    Check every nucleotide mutation that occurred on a branch to determine whether or not it is synonymous.

    For each node, all nucleotide mutations that occurred in parents of the node are applied to the reference sequence to give the genome prior to this node. Then, each nucleotide mutation at the node is made to the appropriate codon from this genome and determined to be synonymous or nonsynonymous.

    Returns a dictionary of synonymous mutations where the key is a gene and the value is a list of synonymous mutations in this gene.
    """
    parent_diffs_pos = [int(k) for k,v in parent_diffs_from_ref.items()]


    # make dictionary of synonymous (and noncoding) mutations to add to tree
    syn_muts = {}

    # don't care about deletions because they are obviously not synonymous
    for mut in nuc_muts_on_branch:
        if mut[-1]!= '-' and mut[0]!='-':
            mut_pos = int(mut[1:-1])
            # find what gene this mut happens in
            if (mut_pos-1) in reference_gene_locations.keys():
                mut_gene = reference_gene_locations[mut_pos-1]
                mut_codon_num = reference_gene_codon[mut_pos-1][0]
                mut_codon_pos = reference_gene_codon[mut_pos-1][1]

                # find the reference sequence of the codon this mutation occurs in
                codon_ref_aa = reference_sequence_aa[mut_gene][mut_codon_num]

                codon_ref_nt = reference_sequence_nt[mut_gene][(mut_codon_num*3):(mut_codon_num*3+3)]

                # check if a mutation occurred within the same codon in a parent
                # and if so, change the reference codon sequence accordingly,
                # to tell whether the mutation at this branch is synonymous or not
                codon_genome_pos = list(range((mut_pos-1-mut_codon_pos),(mut_pos-1-mut_codon_pos+3)))

                parent_codon = codon_ref_nt
                for parent_diff in parent_diffs_pos:
                    parent_diff_zero_based = parent_diff-1
                    if parent_diff_zero_based in codon_genome_pos:
                        parent_diff_pos = codon_genome_pos.index(parent_diff_zero_based)
                        parent_codon = MutableSeq(str(codon_ref_nt))
                        parent_codon[parent_diff_pos] = parent_diffs_from_ref[parent_diff]
                        parent_codon = Seq(parent_codon)


                codon_mutated = MutableSeq(str(parent_codon))
                #if deletion (or seq error) has happened at neighboring nucleotide
                if '-' in codon_mutated:
                    pass
                else:
                    codon_mutated[mut_codon_pos] = mut[-1]
                    codon_mutated = Seq(codon_mutated)
                    codon_mutated_translation = codon_mutated.translate()

                    if str(codon_ref_aa) == str(codon_mutated_translation):
                        if mut_gene in syn_muts.keys():
                            syn_muts[mut_gene] += [mut]
                        else:
                            syn_muts[mut_gene] = [mut]



            else:
                if 'noncoding' in syn_muts.keys():
                    syn_muts['noncoding'] += [mut]
                else:
                    syn_muts['noncoding'] = [mut]

    return syn_muts



def add_syn_mut_attribute(tree, reference_gene_locations, reference_gene_codon, reference_sequence_nt, reference_sequence_aa):
    """
    For each node on the tree, add a node attribute 'syn_muts', which is a dictionary with genes as keys and the value is a list of all synonymous mutations that occur at that node within the gene
    """

    for node in tree.find_clades():

        node.branch_attrs['syn_muts'] = {}

        # only care if this branch has some nucleotide mutations
        if hasattr(node, 'branch_attrs'):
            if 'nuc' in node.branch_attrs['mutations']:

                nuc_muts_on_branch = node.branch_attrs['mutations']['nuc']

                node_path = get_parent(tree, node)

                nucleotide_mut_path = []

                # find all nucleotide mutations that happened in parents,
                # in case they affect codons mutated on this branch
                for parent in node_path[-1]:
                    if hasattr(parent, 'branch_attrs'):
                        if 'nuc' in parent.branch_attrs['mutations']:
                            nucleotide_mut_path+=parent.branch_attrs['mutations']['nuc']

                parent_diffs_from_ref = nuc_changes_from_root(nucleotide_mut_path)

                syn_muts_dict = determine_synonymous(nuc_muts_on_branch, parent_diffs_from_ref, reference_gene_locations, reference_gene_codon, reference_sequence_nt, reference_sequence_aa)

                node.branch_attrs['syn_muts'] = syn_muts_dict
    return tree

def add_changes_from_root_attr(tree, gene_list):
    """
    For each node, find all mutations that occurred between the root and the node. Separate these mutations by gene and by synonymous/nonsynonymous. Add a dictionary of these mutations as an attribute to each node
    """
    for node in tree.find_clades():
        #Find all parents of the node (includes node too)
        parents = get_parent(tree, node)

        #Find mutations that occur in the parents
        parents_nonsyn_muts = {g:[] for g in gene_list}
        parents_syn_muts = {g:[] for g in gene_list}


        for parent in parents:
            if hasattr(parent, "branch_attrs") and "mutations" in parent.branch_attrs:
                for gene in gene_list:
                    if gene in parent.branch_attrs["mutations"]:
                        parents_nonsyn_muts[gene]+=parent.branch_attrs["mutations"][gene]

            if hasattr(parent, 'node_attrs') and 'syn_muts' in parent.branch_attrs:
                for gene in gene_list:
                    if gene in parent.branch_attrs['syn_muts']:
                        parents_syn_muts[gene]+=parent.branch_attrs["syn_muts"][gene]

        # remove reversions from mutation lists
        for g,m in parents_nonsyn_muts.items():
            muts_wo_reversions = remove_reversions(m)
            parents_nonsyn_muts[g] = muts_wo_reversions

        for g,m in parents_syn_muts.items():
            muts_wo_reversions = remove_reversions(m)
            parents_syn_muts[g] = muts_wo_reversions

        # count deletion of adjacent nucleotides as one mutation event
        for g,m in parents_nonsyn_muts.items():
            muts_consolidated_deletions = consolidate_deletions_2(m)
            parents_nonsyn_muts[g] = muts_consolidated_deletions

        for g,m in parents_syn_muts.items():
            muts_consolidated_deletions = consolidate_deletions_2(m)
            parents_syn_muts[g] = muts_consolidated_deletions

        node.branch_attrs["changes_from_root"] = {'nonsyn': parents_nonsyn_muts, 'syn': parents_syn_muts}

    return tree


def remove_reversions(mutation_list):
    """
    If site mutates and then reverts, do not count this in the mutation tally.
    If site mutates and then mutates again (but not a reversion), count only the second mutation
    """
    mutation_list_pos = [int(x[1:-1]) for x in mutation_list]
    sites_mutated_twice = set([x for x in mutation_list_pos if mutation_list_pos.count(x) > 1])
    # find if twice-mutated site was a reversion
    for site in sites_mutated_twice:
        muts_at_site = [mut for mut in mutation_list if int(mut[1:-1])==site]

        # if site reverts, remove all mutations at this site
        if muts_at_site[0][0] == muts_at_site[-1][-1]:
            for mut in range(len(muts_at_site)):
                mutation_list.remove(muts_at_site[mut])
        # if the site mutates multiple times, but doesn't revert, keep last mutation
        else:
            for mut in range(len(muts_at_site)-1):
                mutation_list.remove(muts_at_site[mut])
    return mutation_list


def consolidate_deletions_2(mutation_list):
    """
    Consolidates deletions at adjacent sites into one deletion event. Returns a list of all mutations, including with consolidated deletions
    """

    without_deletions = [x for x in mutation_list if x[-1]!='-' and x[0]!='-']
    #consolidate deletions and reversions
    deletions_only = [x for x in mutation_list if x[-1]=='-' or x[0]=='-']
    deletions_only.sort(key=lambda x:x[1:-1])


    #keep track of start of separate deletions
    separate_deletions = []

    # if there are deletions, count a run of consecutive sites as a single deletion/mutation
    if len(deletions_only) != 0:
        separate_deletions.append(deletions_only[0])

        deletion_tracker = int(deletions_only[0][1:-1])

        for deletion in deletions_only[1:]:

            deleted_pos = int(deletion[1:-1])
            if deleted_pos == deletion_tracker+1:
                pass
            else:
                separate_deletions.append(deletion)
            deletion_tracker = deleted_pos

    consolidated_mutation_list = separate_deletions + without_deletions

    return consolidated_mutation_list


def DateToStr(number):
    """
    Convert decimal date into Mon-Year
    """

    year = int(number)
    d = timedelta(days=(number - year)*365)
    day_one = datetime(year,1,1)
    date = d + day_one

    date_str = date.strftime('%b-%d-%Y')
    return date_str
