# Norovirus
Nextstrain analysis of Norovirus

## Introduction
This is the [Nextstrain build for Norovirus](https://nextstrain.org/community/blab/norovirus/all)

## Data Curation
All sequence data is from Vipr or Genbank. The full Norovirus genomic length is ~7,547 bp long. In this build, we filtered for human Norovirus sequences that are at least 5032bp long (2/3 of the full length). We ended up with a dataset of 1981 sequences from 1968-2022, from 42 countries.
### Obtaining large dataset metadata and sequence output files
1. Break GenomicFastaResults.fasta file (1981 sequences) into 3 using *`seqkit split GenomicFastaResults.fasta -n 703`
      * 703,703, 575 sequences in 3 output files
2. Put 3 files into [norovirus typing tool](https://www.genomedetective.com/app/typingtool/nov/)
3. Concatenate resulting csv metadata files
      * Delete header of files 2 and 3 and concatenate using *`cat file1.csv file2.csv file3.csv > result.csv`
4. Use csvtk to select only relevant columns (ORF2_type)
       *`cat input.csv | csvtk cut -f strain ORF2_type > output.csv`
5. Use regex to delete everything after ‘|’ in strain field
       *`cat input.csv | csvtk replace -f 1 -p "\|.*$" -r "" > output.csv`
6. Convert output file from csv to tsv
7. Join using csvtk
       *`csvtk -t join -f "strain" results/big_metadata.tsv results/genotype_metadata_parsed.tsv --left-join  --na "NA" > metadata_parsed.tsv`
8. Fix Viet_Nam typing issue in output metadata file running fix_vietnam.py script
## Further Reading
Relevant papers for further reading:
* [Norwalk Virus Minor Capsid Protein VP2 Associates within the VP1 Shell Domain](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624303/)
* [Deep Sequencing of Norovirus Genomes Defines Evolutionary Patterns in an Urban Tropical Setting](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4178781/)
