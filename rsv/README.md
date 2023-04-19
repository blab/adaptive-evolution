### Phylogenetic trees of Respiratory Syncytial Virus
This repository builds trees of the Respiratory Syncytial Virus (RSV) subtypes -A and -B from full genome sequences. 
The final trees are viewable at https://nextstrain.org/groups/blab/rsv/A and https://nextstrain.org/groups/blab/rsv/B, and are built in multiple steps using a modified Nextstrain pipeline. 
The steps (outlined below) are used to assign subtype -A or -B to each genome, and to build alignments and trees that are consistent with the large duplication events that occur within the G gene.

#### Step 0: 
- Download all RSV data from ViPR as `vipr_download.fasta`
- Run `format_downloaded_genomes.ipynb` to format dates and output `rsv.fasta`
- Run Snakefile in `rsv_step0/`
- Use `label_rsv_subtypes.ipynb` to separate the unaligned `rsv.fasta` file into RSV-A and RSV-B fasta files based off the tree
- Use `extract_gene_fastas.ipynb` to create separate fasta files for the whole genome, G and F genes

#### Step 1:
- Run Snakefile in `rsv_step1/` to construct separate trees for RSV-A and RSV-B using a reference sequences that contains the G gene duplication

#### Step 2:
- Separate sequences with the duplication from those without using `separate_strains_with_dup.ipynb`
- Run Snakefile in `rsv_step2/` to build two trees for each subtype: one includes only sequences with the G duplication and is aligned to a reference containing the duplication, the other contains only sequences without the G duplication and is aligned to the same reference but with the duplication removed

#### Step 3:
- Use `separate_strains_with_dup.ipynb` to add a column to the `metadata_{subtype}.tsv` files saying whether that strain has the G gene duplication or not (this will be used as a `color_by` option). This also saves the `metadata_{subtype}.tsv` file into `rsv_step3/data/` 
- Run `make_step3_alignment.ipynb` to add placeholders to the aligned sequences that do not have the duplication and to concatenate the alignment files produced in Step 2. This will save `aligned_{subtype}_all.fasta` to `rsv_step3/data/`
- Run `merge_insertion_tsvs.ipynb` to retrieve information about insertions that were inferred during the alignment in Step 2. These will be used to label insertions on branches
- Run Snakefile in `rsv_step3/` to output the final tree
