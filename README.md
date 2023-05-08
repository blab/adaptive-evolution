# An Atlas of Adaptive Evolution in Endemic Human Viruses

**Kathryn Kistler** <sup>1,2</sup>, **Trevor Bedford** <sup>1,2,3</sup> <br />
<sup>1</sup> *Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, United States*<br />
<sup>2</sup> *Howard Hughes Medical Institute, Seattle, WA, United States*

Through antigenic evolution, viruses like influenza A/H3N2 evade recognition by neutralizing antibodies elicited by previous infection or vaccination. This means that a person with antibodies well-tuned to an initial infection will not be protected against the same virus years later and that vaccine-mediated protection will decay. It is not fully understood which of the many endemic human viruses evolve in this fashion and it is the goal of this manuscript to expand that knowledge. To do so, we assess adaptive evolution across the viral genome in 28 endemic viruses, spanning a wide range of viral families and transmission modes. We find that surface proteins consistently show the highest rates of adaptation, and estimate that ten viruses in this panel undergo antigenic evolution to selectively fix mutations that enable the virus to escape recognition by prior immunity. We compare overall rates of amino acid substitution between these antigenically-evolving viruses and SARS-CoV-2, showing that SARS-CoV-2 viruses are accumulating protein-coding changes at substantially faster rates than the endemic viruses. 

---
#### Structure of this repository

This repository contains the data, results, and analysis code to compute rates of adaptation across a panel of human pathogenic viruses. 
>[`adaptive-evolution-analysis/`](https://github.com/blab/adaptive-evolution/tree/master/adaptive-evolution-analysis) : contains the code to estimate rates of adaptative evolution as well as the raw results and Jupyter notebooks that read in these results and generate the figures in the manuscript<br /><br />
>`pathogen_name/`: contains the output of a Nextstrain build that is used an in input for the calculation of rates of adaptation. The necessary files include a FASTA alignment, a metadata file with information about sampling dates, a Genbank reference file, and a Nextstrain tree JSON file.

---
#### Interactive results

The results, highlighting a comparison of of the relative rates of adaptation between viruses as well as within each virus' genome, can be interactively viewed [here](https://blab.github.io/atlas-of-viral-adaptation/)
