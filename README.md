## About

`adaptive-evolution-nextstrain` estimates adaptive evolution in a virus by calculating rates of adaptation using methodology developed in [Bhatt et al 2011](https://pubmed.ncbi.nlm.nih.gov/21415025/) and [Bhatt et al 2010](https://www.sciencedirect.com/science/article/abs/pii/S1567134809001324). A rate of adaptation is calculated for the receptor-binding domain, membrane-fusion domain, and polymerase. Viruses that undergo antigenic drift will exhibit a higher rate of adaptive evolution in the receptor-binding domain than the membrane-fusion domain or polymerase. 

This analysis is designed to be implemented on any virus in coordination with [Nextstrain](https://nextstrain.org), and uses results files of a Nextstrain build as input. 

To run `adaptive-evolution-nextstrain` on a virus:
1. Run Nextstrain build
2. Set up the Config file
3. Run `adaptive-evolution-nextstrain`


## Run Nextstrain build

**If a Nextstrain build already exists for the virus: 
1. Clone the Github repo for the build
2. Run the build using Snakemake

**If a Nextstrain build does not exist for the virus:
1. Set up a new build 
2. Run the build using Snakemake

## Set up the Config file

The config file for each virus should be named `adaptive_evolution_config_VIRUS.json`, and should be located within the `adaptive-evolution/config` directory.

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

## Run `adaptive-evolution-nextstrain`