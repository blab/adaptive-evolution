## About

Adaptive evolution is evolution in response to a selective pressure that increases fitness. In the context of viruses, the rate at which a protein evolves adaptively speaks to its evolutionary potential to undergo antigenic drift (evade antibody recognition), cross-species barriers, or escape drugs. In particular, a high rate of adaptation in a virus that has been endemic in humans for decades indicates that this virus is undergoing continuous adaptation, likely stemming from a changing selective landscape caused by an evolutionary arms race between the virus and its host. Because viral surface proteins are the primary targets of neutralizing antibodies, high rates of adaptation in viral surface protein typically indicates antigenic drift.

This project seeks to compare adaptive evolution across a panel of human pathogenic RNA viruses. Each of these viruses contains a polymerase, which we expect to be relatively conserved, as well as surface proteins (or protein subunits) that bind host cells and fuse host and viral membranes. Using standardized methods, evolution in these proteins can be compared between viruses, with a high rate of adaptive evolution in surface proteins serving as a computational predication of antigenic evolution. 

This project is meant to assess adaptive evolution of a virus using sequence data, and has been built to work in coordination with a Nextstrain pathogen build. This project runs the following analyses on the output of a Nextstrain build  
1. Topology of phylogeny (uses tree)
> Ladder-like trees are indicative of ongoing adaptive evolution, while bushier trees signal a lack of continuous adaptative evolution
2. dN/dS over time (uses tree)
> An excess of nonsynonymous change relative to the expectation indicates selection for functional change. Divergence is calculated as the number of observed mutations (nonsynonymous or synonymous) divided by the number of possible sites (nonsynonymous or synonymous) for every internal branch within a 10-year time window. 
3. Rate of adaptation (uses alignment)
> The rate of adaptation (given in adaptive substitutions per codon per year) can be calculated and used to directly compare viruses. This is calculated directly from the alignment in a sitewise fashion. It defines a neutral class of synonymous and mid-frequency nonsynonymous mutations, and detects adaptive substitutions as an excess of high-frequency or fixed nonsynonymous mutations relative to this neutral class.


To run these analyses on a virus:
1. Run Nextstrain build
2. Set up the Config file
3. Run analyses of adaptive evolution


## Run Nextstrain build

**See [Nextstrain documentation](https://docs.nextstrain.org/) for instructions.

## Set up the Config file

The config file for each virus should be named `adaptive_evo_config_VIRUS.json`, and should be located within the `config/` directory.

The configuration input is written in json format. The following keys are required: 

	"virus": name of virus
	"virus_family": name of virus family
	"subtype": "True" or "False". Whether virus has subtypes. If "True", "subtypes" input must also be specified
	"color": hexcode for color to plot results
	"alignment_file": path to the FASTA alignment file, relative to `bhatt_nextstrain.ipynb`
	"meta_file": path to the TSV metadata file, relative to `bhatt_nextstrain.ipynb`
	"reference_file": path to the Genbank reference genome file, relative to `bhatt_nextstrain.ipynb`
	"ha_protein": dictionary input with the required entry "virus_gene": name of the viral gene equivalent to HA (surface protein with receptor-binding functionality), as specified in the reference Genbank file
	"receptor_binding": dictionary input with the required entry "virus_gene": name of the receptor-binding gene or domain, as specified in the reference Genbank file. If the receptor-binding domain is not listed in the Genbank file, the location of this domain must be specified with "receptor_binding_location" input
	"membrane_fusion": dictionary input with the required entry "virus_gene": name of the membrane-fusion gene or domain, as specified in the reference Genbank file. If the membrane-fusion domain is not listed in the Genbank file, the location of this domain must be specified with "membrane_fusion_location" input
	"polymerase": dictionary input with the required entry "virus_gene": name of the viral polymerase, as specified in the reference Genbank file

The following inputs are optional or conditionally required:
	"subtypes": list of all subtypes. Only required if "subtype": "True"
	"specify_location": use this if one of the required domains (such as "receptor_binding") is a sub-domain of another gene that does not have a separate "alignment_file" or is not listed as a gene in the Genbank "reference_file". an optional key entry in the "ha_protein", "receptor_binding", "membrane_fusion" and "polymerase" dictionaries. The value of "specify_location" is a dictionary with 2 required keys: "parent_gene" and "location". If "subtypes" are specified, then "location_{subtype}" must be supplied for each subtype instead of "location". The "parent_gene" is the gene that contains this domain. The "location" is the the domain's position within the parent gene or the whole genome (whichever format the reference file is in). The format for location entries is "[[start_position, end_position]]" or "[[start_position1, end_position1], [start_position2, end_position2]]" if the domain is non-contiguous.

## Run analyses of adaptive evolution