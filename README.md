# An atlas of continuous adaptive evolution in endemic human viruses

**Kathryn Kistler** <sup>1,2</sup>, **Trevor Bedford** <sup>1,2,3</sup> <br />
<sup>1</sup> *Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, United States*<br />
<sup>2</sup> *Howard Hughes Medical Institute, Seattle, WA, United States*

Through antigenic evolution, viruses such as seasonal influenza evade recognition by neutralizing antibodies. This means that a person with antibodies well tuned to an initial infection will not be protected against the same virus years later and that vaccine-mediated protection will decay. To expand our understanding of which endemic human viruses evolve in this fashion, we assess adaptive evolution across the genome of 28 endemic viruses spanning a wide range of viral families and transmission modes. Surface proteins consistently show the highest rates of adaptation, and ten viruses in this panel are estimated to undergo antigenic evolution to selectively fix mutations that enable the escape of prior immunity. Thus, antibody evasion is not an uncommon evolutionary strategy among human viruses, and monitoring this evolution will inform future vaccine efforts. Additionally, by comparing overall amino acid substitution rates, we show that SARS-CoV-2 is accumulating protein-coding changes at substantially faster rates than endemic viruses.

---
#### Citation

[Kistler KE, Bedford T. 2023. An atlas of continuous adaptive evolution in endemic human viruses. Cell Host Microbe 31: 1-12.](https://doi.org/10.1016/j.chom.2023.09.012)

---

#### Interactive results

The results, highlighting a comparison of of the relative rates of adaptation between viruses as well as within each virus' genome, can be interactively viewed [here](https://blab.github.io/atlas-of-viral-adaptation/)

---

#### Structure of this repository

This repository contains the data, results, and analysis code to compute rates of adaptation across a panel of human pathogenic viruses. 
>[`adaptive-evolution-analysis/`](https://github.com/blab/adaptive-evolution/tree/master/adaptive-evolution-analysis) : contains the code to estimate rates of adaptative evolution as well as the raw results and Jupyter notebooks that read in these results and generate the figures in the manuscript<br /><br />
>`pathogen_name/`: contains the output of a Nextstrain build that is used an in input for the calculation of rates of adaptation. The necessary files include a FASTA alignment, a metadata file with information about sampling dates, a Genbank reference file, and a Nextstrain tree JSON file.
