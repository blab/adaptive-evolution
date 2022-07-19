from __future__ import print_function
from pymol import cmd
import ast


def saveClusters(pdb, virus, subtype, gene, mutation_type, exposure_cutoff, neighbor_cutoff):

    cmd.fetch(pdb)
    #get surface residues (at exposure_cutoff) that have at least one fixation
    if str(subtype)!=str(None):
        fixed_exposed_file = f'/Users/katekistler/nextstrain/adaptive-evolution/adaptive_loci_results/fixations_per_site/results/fixed_exposed_{virus}_{subtype}_{mutation_type}_{exposure_cutoff}.txt'
    else:
        fixed_exposed_file = f'/Users/katekistler/nextstrain/adaptive-evolution/adaptive_loci_results/fixations_per_site/results/fixed_exposed_{virus}_{mutation_type}_{exposure_cutoff}.txt'
    with open(fixed_exposed_file) as f:
        lines = f.readlines()
    surface_residues_w_fixations = []
    for line in lines:
        surface_residues_w_fixations.append(ast.literal_eval(line))

    #keep track of clusters of surface residues with fixations
    clusters = {}

    clusternum_by_residue = {}
    #keep track of all residues already included in a clusters
    already_in_cluster = []

    #assign each separate cluster a number
    cluster_num = 1

    for x in surface_residues_w_fixations:
        chain = x[0]
        residue = x[1]

        #select all surface-exposed residues with fixations
        cmd.select("surface_fixed",f"br. all within {neighbor_cutoff} of (chain {chain} and resi {residue})")

        #find their neighbors
        myspace = {'neighboringresidues': []}
        #use alpha carbon to just get one listing per residue
        cmd.iterate("(surface_fixed) and n. CA", 'neighboringresidues.append((chain, resi))', space=myspace)

        #check if neighbors have fixations
        for neighbor in myspace['neighboringresidues']:
            #self will be called a neighbor, don't want this
            if x!= neighbor:
                if neighbor in surface_residues_w_fixations:
                    #find whether either x or neighbor is already part of a cluster to add on to
                    if x in already_in_cluster:
                        this_cluster = clusternum_by_residue[x]
                        already_in_cluster.append(neighbor)
                    elif neighbor in already_in_cluster:
                        this_cluster = clusternum_by_residue[neighbor]
                        already_in_cluster.append(x)
                    #otherwise, start a new cluster
                    else:
                        this_cluster = cluster_num
                        cluster_num+=1
                        already_in_cluster.append(x)
                        already_in_cluster.append(neighbor)

                    #store information of what cluster each residue is in
                    if x not in clusternum_by_residue.keys():
                        clusternum_by_residue[x] = this_cluster
                    if neighbor not in clusternum_by_residue.keys():
                        clusternum_by_residue[neighbor] = this_cluster

                    if this_cluster in clusters.keys():
                        if x not in clusters[this_cluster]:
                            clusters[this_cluster].append(x)
                        if neighbor not in clusters[this_cluster]:
                            clusters[this_cluster].append(neighbor)
                    else:
                        clusters[this_cluster] = [x, neighbor]



    if str(subtype)!=str(None):
        output_filename = f'{virus}_{subtype}_{mutation_type}_{exposure_cutoff}_{neighbor_cutoff}_clusters'
    else:
        output_filename = f'{virus}_{mutation_type}_{exposure_cutoff}_{neighbor_cutoff}_clusters'
    with open(f"/Users/katekistler/nextstrain/adaptive-evolution/adaptive_loci_results/fixations_per_site/results/clusters/{output_filename}.txt", "w") as f:
        f.write(str(clusters))
        f.write('\n')

    f.close()


cmd.extend("saveClusters", saveClusters)
