import math
import json
import random
import ast
import re
from os import path
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

def readin_virus_config(virus):
    config_json = f'config/adaptive_evo_config_{virus}.json'
    with open(config_json) as json_handle:
        configs = json.load(json_handle)

    return configs


def frequency_binning(x, midfreq_high, midfreq_low):
    """
    given a polymorphism frequency, return bin
    """
    #nan frequencies are when there is no sequence coverage at the given position
    if math.isnan(x):
        f_bin = float('nan')
    else:
        if x == 1.0:
            f_bin = 'f'
        elif x>=midfreq_high:
            f_bin = 'h'
        elif x<midfreq_high and x>=midfreq_low:
            f_bin = 'm'
        elif x<midfreq_low:
            f_bin='l'

    return f_bin


def walk_through_sites(outgroup_seq, outgroup_aa_seq, alignment_seqs, midfreq_high, midfreq_low):

    #at each site, count number of viruses with polymorphism
    count_polymorphic = np.zeros(len(outgroup_seq))
    #at each site, count totaly number of viruses
    count_total_unambiguous = np.zeros(len(outgroup_seq))

    count_replacement_mutations = np.zeros(len(outgroup_seq))
    count_silent_mutations = np.zeros(len(outgroup_seq))

    #at each site, list of nucleotide from each virus
    ingroup_bases = [[] for x in range(len(outgroup_seq))]

    for seq in alignment_seqs:
        if len(seq) != len(outgroup_seq):
            print(len(seq), len(outgroup_seq))
            print(seq)
        elif len(seq) == len(outgroup_seq):
            for pos in range(len(outgroup_seq)):
                outgroup_nt = str(outgroup_seq[pos])
                virus_nt = str(seq[pos])

                #skip ambiguous sites
#                 if virus_nt != 'N' and virus_nt != 'W':
#                     if outgroup_nt != 'N' and outgroup_nt != 'W':
                if virus_nt in ['A', 'C', 'G', 'T']:
                    if outgroup_nt in ['A', 'C', 'G', 'T']:
                        ingroup_bases[pos].append(virus_nt)
                        count_total_unambiguous[pos]+=1
                        if virus_nt != outgroup_nt:
                            count_polymorphic[pos]+=1
                            #determine silent or replacement
                            codon = math.floor(pos/3)

                            codon_pos = pos-(codon*3)
                            if codon_pos == 0:
                                codon_nt = virus_nt+outgroup_seq[pos+1:(pos+3)]
                            elif codon_pos == 1:
                                codon_nt = outgroup_seq[pos-1]+virus_nt+outgroup_seq[pos+1]
                            elif codon_pos == 2:
                                codon_nt = outgroup_seq[(pos-2):(pos)]+virus_nt


                            if isinstance(codon_nt, str):
                                codon_nt = Seq(codon_nt)

                            codon_aa = codon_nt.translate()


                            outgroup_aa = outgroup_aa_seq[codon]

                            if outgroup_aa != 'X':
                                if codon_aa != outgroup_aa:
                                    count_replacement_mutations[pos]+=1
                                elif codon_aa == outgroup_aa:
                                    count_silent_mutations[pos]+=1

    polymorphic_frequencies = count_polymorphic/count_total_unambiguous

    replacement_score = count_replacement_mutations/count_total_unambiguous

    freq_bins = [frequency_binning(x, midfreq_high, midfreq_low) for x in polymorphic_frequencies]

    return freq_bins, replacement_score, ingroup_bases

def determine_site_type(outgroup, ingroup):

    ingroup_bases_nan = set(ingroup)
    #remove 'nan's
    ingroup_bases = {x for x in ingroup_bases_nan if pd.notna(x)}


    if len(ingroup_bases) == 0:
        site_type = None

    elif len(ingroup_bases) != 0:
        #all ingroup bases are identical
        if len(ingroup_bases) == 1:
            if outgroup in ingroup_bases:
                site_type = 1
            elif outgroup not in ingroup_bases:
                site_type = 2

        #2 different bases in ingroup
        elif len(ingroup_bases) == 2:
            if outgroup in ingroup_bases:
                site_type = 3
            elif outgroup not in ingroup_bases:
                site_type = 4

        #3 different bases in ingroup
        elif len(ingroup_bases) == 3:
            if outgroup in ingroup_bases:
                site_type = 5
            elif outgroup not in ingroup_bases:
                site_type = 6

        #4 different bases in ingroup
        elif len(ingroup_bases) == 4:
            site_type = 7


    return site_type


def fixation_polymorphism_score(outgroup, ingroup):
    site_type = determine_site_type(outgroup, ingroup)


    if site_type == None:
        Fi = float('nan')
        Pi = float('nan')
    if site_type == 1:
        Fi = 0
        Pi = 0
    elif site_type == 2:
        Fi = 1
        Pi = 0
    elif site_type in [3,5,7]:
        Fi = 0
        Pi = 1
    elif site_type == 4:
        Fi = 0.5
        Pi = 0.5
    elif site_type == 6:
        Fi = (1/3)
        Pi = (2/3)

    return Fi, Pi


def assign_fi_pi(outgroup_seq, ingroup_bases):

    #at each site, record Fi
    Fi_all = np.zeros(len(outgroup_seq))

    #at each site, record Pi
    Pi_all = np.zeros(len(outgroup_seq))

    for pos in range(len(outgroup_seq)):
        outgroup_nt = outgroup_seq[pos]
        ingroup_nts = ingroup_bases[pos]
        Fi, Pi = fixation_polymorphism_score(outgroup_nt, ingroup_nts)
        Fi_all[pos] = Fi
        Pi_all[pos] = Pi

    return Fi_all,


def subset_viruses_nextstrain_build(virus, subtype, gene, window, min_seqs, year_max, year_min):

    configs = readin_virus_config(virus)
    standard_gene = standardize_gene_name(virus, gene)


    #Find reference, alignment and meta files (some sub-genic regions may use files from a gene or a whole genome)
    if 'specify_location' in configs[standard_gene].keys():
        parent_gene = configs[standard_gene]['specify_location']['parent_gene']
        reference_file = configs['reference_file'].format(virus=virus, subtype=subtype, gene=parent_gene)
        alignment_file = configs['alignment_file'].format(virus=virus, subtype=subtype, gene=parent_gene)
        meta_file = configs['meta_file'].format(virus=virus, subtype=subtype, gene=parent_gene)
        #some are comma-separated, some are tab-separated
        metafile_sep = configs['metafile_sep']
    else:
        reference_file = configs['reference_file'].format(virus=virus, subtype=subtype, gene=gene)
        alignment_file = configs['alignment_file'].format(virus=virus, subtype=subtype, gene=gene)
        meta_file = configs['meta_file'].format(virus=virus, subtype=subtype, gene=gene)
        metafile_sep = configs['metafile_sep']

    #Find gene location, if domain is sub-genic or reference file contains multiple genes
    gene_location = False
    #If domain is sub-genic, fetch its position (within genome or parent gene) from config file


    if 'specify_location' in configs[standard_gene].keys():
        if subtype==None:
            gene_location_key = "location"
        else:
            gene_location_key = "location_"+str(subtype)

        gene_location_list = ast.literal_eval(configs[standard_gene]['specify_location'][gene_location_key])
        #Need to deal with domains the are not contiguous
        if len(gene_location_list)==1:
            gene_location = SeqFeature(FeatureLocation(gene_location_list[0][0], gene_location_list[0][1]))
        else:
            compound_locations = []
            for location in gene_location_list:
                compound_locations.append(FeatureLocation(location[0], location[1]))
            gene_location = CompoundLocation(compound_locations)

    #Find gene location from reference files
    else:
        for seq_record in SeqIO.parse(reference_file, "genbank"):
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    if 'gene' in feature.qualifiers.keys():
                        if feature.qualifiers['gene'][0].lower() == gene.lower():
                            gene_location = feature.location
                    elif feature.qualifiers['product'][0].lower() == gene.lower():
                        gene_location = feature.location

    #Subset data based on time windows
    meta = pd.read_csv(meta_file, sep = metafile_sep)
    meta.drop(meta[meta['date']=='?'].index, inplace=True)
    meta.dropna(subset=['date'], inplace=True)
    meta['year'] = meta['date'].str[:4].astype('int')
    if year_max:
        meta.drop(meta[meta['year']>year_max].index, inplace=True)
    if year_min:
        meta.drop(meta[meta['year']<year_min].index, inplace=True)

    date_range = meta['year'].max() - meta['year'].min()
    #Remove egg- and cell-passaged strains
    meta.drop(meta[meta['strain'].str[-4:]=='-egg'].index, inplace=True)
    meta.drop(meta[meta['strain'].str[-5:]=='-cell'].index, inplace=True)

    #Limit meta data to only strains in alignment file
    aligned_isolates = []
    with open(alignment_file, "r") as aligned_handle:
        for isolate in SeqIO.parse(aligned_handle, "fasta"):
            aligned_isolates.append(isolate.id)
    aligned_isolates_df = pd.DataFrame(aligned_isolates, columns=['strain'])
    meta = meta.merge(aligned_isolates_df, on='strain', how='inner')


    #Group viruses by time windows
    virus_time_subset = {}
    if window == 'all':
        years = str(meta['year'].min()) + '-' + str(meta['year'].max())
        virus_time_subset[years] = meta['strain'].tolist()
    else:
        date_window_start = meta['year'].min()
        date_window_end = meta['year'].min() + window
        while date_window_end <= meta['year'].max():
            years = str(date_window_start) + '-' + str(date_window_end)
            strains = meta[(meta['year']>=date_window_start) & (meta['year']<date_window_end)]['strain'].tolist()
            virus_time_subset[years] = strains


            #sliding window
            date_window_end += 1
            date_window_start += 1


    #Only use time points with enough data:
    virus_time_subset = {k:v for k,v in virus_time_subset.items() if len(v)>=min_seqs}

    year_windows = []
    seqs_in_window = []

    #Find outgroup sequence from strains at first time point(to make consensus from)
    first_window = True
    first_window_strains = []
    first_window_sequences = []

    alignment_time_subset = {}


    for years, subset_viruses in virus_time_subset.items():

        year_windows.append(years)
        seqs_in_window.append(len(subset_viruses))
        alignment_time_subset[years] = []

        #make consensus sequence at first time point
        if first_window == True:
            first_window_strains+=subset_viruses
            first_window = False


        with open(alignment_file, "r") as aligned_handle:
            for isolate in SeqIO.parse(aligned_handle, "fasta"):
                if isolate.id in first_window_strains:
                    if gene_location:
                        gene_record = SeqRecord(seq = gene_location.extract(isolate.seq),
                                                id = isolate.id, description = gene)
                    else:
                        gene_record = SeqRecord(seq = isolate.seq,
                                                id = isolate.id, description = gene)
                    first_window_sequences.append(gene_record)
                if isolate.id in subset_viruses:
                    if gene_location:
                        alignment_time_subset[years].append(gene_location.extract(isolate.seq))
                    else:
                        alignment_time_subset[years].append(isolate.seq)

    first_window_alignment = MultipleSeqAlignment(first_window_sequences)
    outgroup_seq = AlignInfo.SummaryInfo(first_window_alignment).gap_consensus(ambiguous ='N')
    outgroup_seq_aa = outgroup_seq.translate()

    return virus_time_subset, alignment_time_subset, outgroup_seq, outgroup_seq_aa, year_windows, seqs_in_window



def bootstrap_alignment(bootstrap_codon_order, sequences):
    """
    for each time point, create sample alignment of same size as emperical alignment
    """

    bootstrap_alignment_seqs = []
    for virus_seq in sequences:
        virus_seq_str = str(virus_seq)
        virus_codons = [virus_seq_str[i:i+3] for i in range(0, len(virus_seq_str), 3)]
        bootstrap_virus = ''.join([virus_codons[x] for x in bootstrap_codon_order])
        bootstrap_alignment_seqs.append(bootstrap_virus)

    return bootstrap_alignment_seqs


def bootstrap_ancestral(outgroup_seq):
    """
    sample codons from emperical ancestral sequence with replacement
    """
    outgroup_seq_str = str(outgroup_seq)
    #sample codons with replacement
    ancestral_codons = [outgroup_seq_str[i:i+3] for i in range(0, len(outgroup_seq_str), 3)]
    bootstrap_codon_order = random.choices(range(len(ancestral_codons)), k=len(ancestral_codons))
    bootstrap_ancestral_seq = ''.join([ancestral_codons[x] for x in bootstrap_codon_order])
    bootstrap_ancestral_seq = Seq(bootstrap_ancestral_seq)
    return bootstrap_ancestral_seq, bootstrap_codon_order

def make_bootstrap_dataset(outgroup_seq, alignment_time_subset):

    bootstrap_ancestral_seq, bootstrap_codon_order = bootstrap_ancestral(outgroup_seq)
    bootstrap_ancestral_seq_aa = bootstrap_ancestral_seq.translate()

    bootstrap_alignment_seqs = {}
    for years, sequences in alignment_time_subset.items():
        bootstrap_sequences = bootstrap_alignment(bootstrap_codon_order, sequences)
        bootstrap_alignment_seqs[years] = bootstrap_sequences


    return bootstrap_ancestral_seq, bootstrap_ancestral_seq_aa, bootstrap_alignment_seqs


def calc_site_stats(alignment_sequences, outgroup_seq, outgroup_aa_seq, midfreq_high, midfreq_low):

    #Find percent polymorphism at each site
    #Also determine whether polymorphism is silent or replacement


    #initiate lists to record all time windows
    frequency_bins = []
    fixation_scores = []
    polymorphism_scores = []
    replacement_scores = []
    silent_scores = []


    for years, alignment_seqs in alignment_sequences.items():

        #calculate stats for each window separately
        freq_bins, replacement_score, ingroup_bases = walk_through_sites(outgroup_seq, outgroup_aa_seq,
                                                                         alignment_seqs,
                                                                         midfreq_high, midfreq_low)
        Fi_all, Pi_all = assign_fi_pi(outgroup_seq, ingroup_bases)
        silent_score = 1-replacement_score

        frequency_bins.append(freq_bins)
        fixation_scores.append(Fi_all)
        polymorphism_scores.append(Pi_all)
        replacement_scores.append(replacement_score)
        silent_scores.append(silent_score)



    return frequency_bins, fixation_scores, polymorphism_scores, replacement_scores, silent_scores


def calc_m_ratio(virus, subtype, gene, window, min_seqs, midfreq_high, midfreq_low, bootstrap, year_max, year_min):
    """
    M=rm/sm
    not expected to vary through time provided that long-term effective population sizes remain sufficiently large
    For each gene, calculate M by combining site count among time points
    """

    configs = readin_virus_config(virus)
    nonantigenic_gene = configs['membrane_fusion']['virus_gene']


    if standardize_gene_name(virus, gene) =='ha_protein' or standardize_gene_name(virus, gene) =='receptor_binding':
        (virus_time_subset, alignment_time_subset,
         outgroup_seq, outgroup_aa_seq,
         year_windows, seqs_in_window) = subset_viruses_nextstrain_build(virus, subtype, nonantigenic_gene, 'all',
                                                        min_seqs, year_max, year_min)
        if bootstrap:
            (bootstrap_ancestral_seq, bootstrap_ancestral_seq_aa,
             bootstrap_alignment_seqs) = make_bootstrap_dataset(outgroup_seq, alignment_time_subset)


    else:
        (virus_time_subset, alignment_time_subset,
         outgroup_seq, outgroup_aa_seq,
         year_windows, seqs_in_window) = subset_viruses_nextstrain_build(virus, subtype, gene, 'all',
                                                        min_seqs, year_max, year_min)
        if bootstrap:
            (bootstrap_ancestral_seq, bootstrap_ancestral_seq_aa,
             bootstrap_alignment_seqs) = make_bootstrap_dataset(outgroup_seq, alignment_time_subset)

    if bootstrap:
            (frequency_bins,
             fixation_scores, polymorphism_scores,
             replacement_scores, silent_scores) = calc_site_stats(bootstrap_alignment_seqs,
                                                                  bootstrap_ancestral_seq, bootstrap_ancestral_seq_aa,
                                                                  midfreq_high, midfreq_low)
    else:
        (frequency_bins,
         fixation_scores, polymorphism_scores,
         replacement_scores, silent_scores) = calc_site_stats(alignment_time_subset,
                                                              outgroup_seq, outgroup_aa_seq, midfreq_high, midfreq_low)


    sm = 0
    rm = 0

    for site in range(len(frequency_bins[0])):
        freq_bin = frequency_bins[0][site]
        if freq_bin == 'm':
            sm+= (polymorphism_scores[0][site]*silent_scores[0][site])
            rm+= (polymorphism_scores[0][site]*replacement_scores[0][site])

    if sm ==0:
        sm = 0.00000000000000001
    m_ratio = rm/sm

    return m_ratio


def bhatt_estimators(gene, outgroup_seq, frequency_bins, year_windows, fixation_scores, polymorphism_scores, replacement_scores, silent_scores, m_ratio):


    #Initiate lists to store a values
    window_midpoint = []
    adaptive_substitutions = []

    #for each window, calculate bhatt estimators
    for years_window in range(len(frequency_bins)):
        window_start = int(year_windows[years_window][0:4])
        window_end = int(year_windows[years_window][-4:])
        window_midpoint.append((window_start + window_end)/2)

        sf = 0
        rf = 0
        sh = 0
        rh = 0
        sm = 0
        rm = 0
        sl = 0
        rl = 0

        #calculate number of sites in different catagories (defined by polymorphic freq at that site)
        window_freq_bins = frequency_bins[years_window]
        for site in range(len(window_freq_bins)):
            freq_bin = window_freq_bins[site]
            #ignore sites with no polymorphisms?
            if freq_bin!='nan':
                if freq_bin == 'f':
                    sf+= (fixation_scores[years_window][site]*silent_scores[years_window][site])
                    rf+= (fixation_scores[years_window][site]*replacement_scores[years_window][site])
                elif freq_bin == 'h':
                    sh+= (polymorphism_scores[years_window][site]*silent_scores[years_window][site])
                    rh+= (polymorphism_scores[years_window][site]*replacement_scores[years_window][site])
                elif freq_bin == 'm':
                    sm+= (polymorphism_scores[years_window][site]*silent_scores[years_window][site])
                    rm+= (polymorphism_scores[years_window][site]*replacement_scores[years_window][site])
                elif freq_bin == 'l':
                    sl+= (polymorphism_scores[years_window][site]*silent_scores[years_window][site])
                    rl+= (polymorphism_scores[years_window][site]*replacement_scores[years_window][site])


#         print(year_windows[years_window])
#         print(sf, rf, sh, rh, sm, rm, sl, rl)

        #Calculate equation 1: number of nonneutral sites
        al = rl - sl*m_ratio
        ah = rh - sh*m_ratio
        af = rf - sf*m_ratio

        #set negative a values to zero
        if al < 0:
            al = 0
        if ah < 0:
            ah = 0
        if af < 0:
            af = 0

#             print(al, ah, af)

        #Calculate the number and proportion of all fixed or high-freq sites that have undergone adaptive change
        number_adaptive_substitutions = af + ah
        adaptive_substitutions.append(number_adaptive_substitutions)
#         proportion_adaptive_sites = (af + ah)/(rf +rh)


    gene_length = len(outgroup_seq)
    adaptive_substitutions_per_codon = [x/gene_length for x in adaptive_substitutions]

    if len(window_midpoint)!=0:
        rate_of_adaptation, intercept, r_value, p_value, std_err = stats.linregress(window_midpoint, adaptive_substitutions_per_codon)
    else:
        rate_of_adaptation = 0

    return window_midpoint, adaptive_substitutions, adaptive_substitutions_per_codon, rate_of_adaptation


def calc_bhatt_a(virus, subtype, gene, window, min_seqs, midfreq_high, midfreq_low, bootstrap, year_max, year_min):
    #Get virus subset
    (virus_time_subset, alignment_time_subset,
     outgroup_seq, outgroup_aa_seq, year_windows, seqs_in_window) = subset_viruses_nextstrain_build(virus, subtype, gene,
                                                                                   window, min_seqs,
                                                                                   year_max, year_min)
#     print(alignment_time_subset, [len(alignment_time_subset[x]) for x in alignment_time_subset.keys()], seqs_in_window)

    #calculate m ratio
    m_ratio = calc_m_ratio(virus, subtype, gene, window,
                           min_seqs, midfreq_high, midfreq_low, False,
                           year_max, year_min)

    #Calculate frequencies for emperical data
    (frequency_bins,
     fixation_scores, polymorphism_scores,
     replacement_scores, silent_scores) = calc_site_stats(alignment_time_subset, outgroup_seq,
                                                          outgroup_aa_seq, midfreq_high, midfreq_low)



    #calculate bhatt estimators
    (window_midpoint, adaptive_substitutions,
     adaptive_substitutions_per_codon, rate_of_adaptation) = bhatt_estimators(gene, outgroup_seq,
                                                                              frequency_bins, year_windows,
                                                                              fixation_scores, polymorphism_scores,
                                                                              replacement_scores, silent_scores, m_ratio)


    n_bootstraps = 100
    bootstrap_count = 0

    bootstrap_adaptive_substitutions = []
    bootstrap_adaptive_substitutions_per_codon = []
    bootstrap_rate_of_adaptation = []
    if bootstrap:
        while bootstrap_count < n_bootstraps:
            bootstrap_count+=1
            #Get bootstrapped ancestral seq and alignment
            (bootstrap_ancestral_seq, bootstrap_ancestral_seq_aa,
             bootstrap_alignment_seqs) = make_bootstrap_dataset(outgroup_seq, alignment_time_subset)


            #Calculate frequencies for bootstrap data
            (bootstrap_frequency_bins,
             bootstrap_fixation_scores, bootstrap_polymorphism_scores,
             bootstrap_replacement_scores, bootstrap_silent_scores) = calc_site_stats(bootstrap_alignment_seqs,
                                                                                      bootstrap_ancestral_seq,
                                                                                      bootstrap_ancestral_seq_aa,
                                                                                      midfreq_high, midfreq_low)
            #Calculate m ratio
            bootstrap_m_ratio = calc_m_ratio(virus, subtype, gene, window,
                                             min_seqs, midfreq_high, midfreq_low, True,
                                             year_max, year_min)

            #calculate bhatt estimators
            (bs_window_midpoint, bs_adaptive_substitutions,
             bs_adaptive_substitutions_per_codon,
             bs_rate_of_adaptation) = bhatt_estimators(gene, bootstrap_ancestral_seq,
                                                              bootstrap_frequency_bins, year_windows,
                                                              bootstrap_fixation_scores,
                                                              bootstrap_polymorphism_scores,
                                                              bootstrap_replacement_scores, bootstrap_silent_scores,
                                                              bootstrap_m_ratio)
            #add these bootstrap values to list
            bootstrap_adaptive_substitutions.append(bs_adaptive_substitutions)
            bootstrap_adaptive_substitutions_per_codon.append(bs_adaptive_substitutions_per_codon)
            bootstrap_rate_of_adaptation.append(bs_rate_of_adaptation)

    if bootstrap:
        return window_midpoint, adaptive_substitutions, adaptive_substitutions_per_codon, rate_of_adaptation, bootstrap_adaptive_substitutions, bootstrap_adaptive_substitutions_per_codon, bootstrap_rate_of_adaptation

    else:
        return window_midpoint, adaptive_substitutions, adaptive_substitutions_per_codon, rate_of_adaptation


def standardize_gene_name_reverse(virus, gene):

    configs = readin_virus_config(virus)

    genes = ['polymerase', 'receptor_binding', 'membrane_fusion']
    gene_names = {x: configs[x]['virus_gene'] for x in genes}

    return gene_names[gene]


def standardize_gene_name(virus, gene):
    configs = readin_virus_config(virus)

    genes = ['polymerase', 'receptor_binding', 'membrane_fusion']
    gene_names = {configs[x]['virus_gene']:x for x in genes}

    return gene_names[gene]




def get_data_to_plot(virus, subtype, gene, bootstrap, window, min_seqs, midfreq_high, midfreq_low, year_max, year_min):

    data_to_plot = []

    if subtype==None:
        virus_subtype = virus
        virus_and_subtype = virus
    else:
        virus_subtype = subtype
        virus_and_subtype = virus+'_'+subtype

    if bootstrap:
        save_json_name = 'bhatt_results_nextstrain/'+str(virus_and_subtype)+'_'+str(gene)+'_bhatt_analysis_bootstrapped.json'
        if path.exists(save_json_name):
            with open(save_json_name) as json_handle:
                json_dict = json.load(json_handle)
                (window_midpoint, adaptive_substitutions,
                 adaptive_substitutions_per_codon,
                 rate_of_adaptation, bootstrap_adaptive_substitutions,
                 bootstrap_adaptive_substitutions_per_codon,
                 bootstrap_rate_of_adaptation) = (json_dict['window_midpoint'],
                                                  json_dict['adaptive_substitutions'],
                                                  json_dict['adaptive_substitutions_per_codon'],
                                                  json_dict['rate_of_adaptation'],
                                                  json_dict['bootstrap_adaptive_substitutions'],
                                                  json_dict['bootstrap_adaptive_substitutions_per_codon'],
                                                  json_dict['bootstrap_rate_of_adaptation'])

        else:

            (window_midpoint, adaptive_substitutions,
             adaptive_substitutions_per_codon,
             rate_of_adaptation, bootstrap_adaptive_substitutions,
             bootstrap_adaptive_substitutions_per_codon,
             bootstrap_rate_of_adaptation) = calc_bhatt_a(virus, subtype, gene, window,
                                                          min_seqs, midfreq_high,
                                                          midfreq_low, bootstrap, year_max, year_min)

            save_json = {'virus': virus, 'subtype':subtype, 'gene': gene, 'window':window, 'min_seqs': min_seqs,
                         'midfreq_high': midfreq_high, 'midfreq_low': midfreq_low,
                         'window_midpoint':window_midpoint, 'adaptive_substitutions':adaptive_substitutions,
                         'adaptive_substitutions_per_codon':adaptive_substitutions_per_codon, 'rate_of_adaptation': rate_of_adaptation,
                         'bootstrap_adaptive_substitutions': bootstrap_adaptive_substitutions,
                         'bootstrap_adaptive_substitutions_per_codon': bootstrap_adaptive_substitutions_per_codon,
                         'bootstrap_rate_of_adaptation':bootstrap_rate_of_adaptation}
            with open(save_json_name, 'w') as outfile:
                json.dump(save_json, outfile)

        slope_sci = rate_of_adaptation * (10**3)
        bs_slope_sci = [x * (10**3) for x in bootstrap_rate_of_adaptation]
        lower_95ci = np.percentile(sorted(bs_slope_sci), 2.5)
        upper_95ci = np.percentile(sorted(bs_slope_sci), 97.5)

        data_to_plot.append({'virus': virus, 'subtype': subtype, 'virus_and_subtype': virus_and_subtype,
                             'gene': standardize_gene_name(virus, gene),
                             'adaptive_subs_per_codon_per_year': slope_sci,
                             'lower_95ci': lower_95ci, 'upper_95ci': upper_95ci,
                             'ci': [lower_95ci, upper_95ci]})




    else:
        save_json_name = 'bhatt_results_nextstrain/'+str(virus_and_subtype)+'_'+str(gene)+'_bhatt_analysis.json'
        if path.exists(save_json_name):
            with open(save_json_name) as json_handle:
                json_dict = json.load(json_handle)
                (window_midpoint, adaptive_substitutions,
                 adaptive_substitutions_per_codon,
                 rate_of_adaptation) = (json_dict['window_midpoint'],
                                        json_dict['adaptive_substitutions'],
                                        json_dict['adaptive_substitutions_per_codon'],
                                        json_dict['rate_of_adaptation'])


        else:
            (window_midpoint, adaptive_substitutions,
             adaptive_substitutions_per_codon,
             rate_of_adaptation) = calc_bhatt_a(virus, subtype, gene, window, min_seqs,
                                                midfreq_high, midfreq_low,
                                                bootstrap, year_max, year_min)


            save_json = {'virus': virus, 'subtype':subtype, 'gene': gene, 'window':window, 'min_seqs': min_seqs,
                         'midfreq_high': midfreq_high, 'midfreq_low': midfreq_low,
                         'window_midpoint':window_midpoint, 'adaptive_substitutions':adaptive_substitutions,
                         'adaptive_substitutions_per_codon':adaptive_substitutions_per_codon,
                         'rate_of_adaptation': rate_of_adaptation}
            with open(save_json_name, 'w') as outfile:
                json.dump(save_json, outfile)

        slope_sci = rate_of_adaptation * (10**3)
        data_to_plot.append({'virus': virus, 'subtype': subtype, 'virus_and_subtype':virus_and_subtype,
                             'gene': standardize_gene_name(virus, gene),
                             'adaptive_subs_per_codon_per_year': slope_sci})


    return data_to_plot


def compare_viruses_adaptive_rate(viruses, standard_genes=['polymerase', 'receptor_binding', 'membrane_fusion'],
                                  window=3, min_seqs=2, bootstrap=False,
                                  midfreq_high=0.75, midfreq_low=0.15, year_max=None, year_min=None, filename=None):



    data_to_plot = []
    color_map = {}
    virus_families = {}

    for virus in viruses:

        configs = readin_virus_config(virus)
        genes = [standardize_gene_name_reverse(virus, x) for x in standard_genes]

        if configs['subtype']=='True':
            subtypes = configs['subtypes']
            for subtype in subtypes:
                virus_and_sub = virus+'_'+subtype
                color_map[virus_and_sub] = configs['color'][subtype]
                virus_families[virus_and_sub] = configs['virus_family']
                for gene in genes:
                    data_to_plot+=get_data_to_plot(virus, subtype, gene, bootstrap, window,
                                                         min_seqs, midfreq_high, midfreq_low, year_max, year_min)

        else:
            subtype=None
            color_map[virus] = configs['color']
            virus_families[virus] = configs['virus_family']
            for gene in genes:
                data_to_plot+=get_data_to_plot(virus, subtype, gene, bootstrap, window,
                                                     min_seqs, midfreq_high, midfreq_low, year_max, year_min)


    df_to_plot = pd.DataFrame(data_to_plot)

    viruses_and_subtypes = list(df_to_plot['virus_and_subtype'].unique())

    x_coords = {}

    all_x_ticks = []
    last_coord = 0.0
    spacing = {**{x:0.25 for x in range(1,8)}, **{y:0.4 for y in range(8,12)}, **{z:0.7 for z in range(12,20)}}
    spacing_genes = {**{x:1.0 for x in range(1,8)}, **{y:1.8 for y in range(8,12)}, **{z:2.4 for z in range(12,20)}}
    for gene in standard_genes:
        x_coords[gene] = {}
        for virus_subtype in viruses_and_subtypes:
            last_coord+=spacing[len(viruses_and_subtypes)]
            x_coords[gene][virus_subtype] = last_coord
            all_x_ticks.append(last_coord)
        last_coord+=spacing_genes[len(viruses_and_subtypes)]

    fig, ax = plt.subplots(figsize=(12,8))

    x_labels = []
    gene_ticks = []


    for gene in standard_genes:
        gene_coords = list(x_coords[gene].values())
        gene_ticks.append(sum(gene_coords)/len(gene_coords))
        x_labels.append(gene)
        for virus_subtype in viruses_and_subtypes:
            x = x_coords[gene][virus_subtype]
            df_row = df_to_plot[(df_to_plot['gene']==gene)&(df_to_plot['virus_and_subtype']==virus_subtype)]
            y = float(df_row['adaptive_subs_per_codon_per_year'])
            if bootstrap:
                err_lower = float(df_row['lower_95ci'])
                err_upper = float(df_row['upper_95ci'])
                ax.vlines( x, err_lower, err_upper, color='black')
            ax.plot(x, y, 'o', ms=12, color=color_map[virus_subtype])




    plt.xticks(gene_ticks, x_labels)

    legend_markers = []
    used_families = []

    for virus_subtype in viruses_and_subtypes:
        if virus_families[virus_subtype] not in used_families:
            legend_markers.append(mlines.Line2D([0], [0], color='w', markerfacecolor='w', marker='o',
                            markersize=12, label=virus_families[virus_subtype]))
            used_families.append(virus_families[virus_subtype])
        legend_markers.append(mlines.Line2D([0], [0], color='w', markerfacecolor=color_map[virus_subtype], marker='o',
                                            markersize=12, label=virus_subtype))
    plt.legend(handles=legend_markers, loc='upper right')

    plt.ylabel('adaptive subs per codon per year (x10^-3)')
#     plt.xlabel('gene')

    # remove box around plot
    sns.despine(left=False, bottom=False)




    if filename:
        fig.savefig(filename, dpi=300)
