## About

Adaptive evolution is evolution in response to a selective pressure that increases fitness. In the context of viruses, the rate at which a protein evolves adaptively speaks to its evolutionary potential to undergo antigenic drift (evade antibody recognition), cross-species barriers, or escape drugs. In particular, a high rate of adaptation in a virus that has been endemic in humans for decades indicates that this virus is undergoing continuous adaptation, likely stemming from a changing selective landscape caused by an evolutionary arms race between the virus and its host. Because viral surface proteins are the primary targets of neutralizing antibodies, high rates of adaptation in viral surface protein typically indicates antigenic drift.

This project seeks to compare adaptive evolution across a panel of human pathogenic RNA viruses. All of these viruses contain a surface proteins (or protein subunit) that binds host cells, and thus is a prime target for the adaptive immune system to block infection. Using standardized methods, evolution in these proteins can be compared between viruses, with a high rate of adaptive evolution in receptor-binding protein serving as a computational predication of antigenic evolution. 

This project is meant to assess adaptive evolution of a virus using sequence data, and has been built to work in coordination with a Nextstrain pathogen build. The analysis code here has been written to consider the evolution of a single lineage of virus that has been sequenced over time, and calculates a rate of adaptation for a given viral protein using a McDonald-Kreitman-based method. This analysis can be done either from a constant outgroup (determined from the consensus sequence at the first timepoint), or from an updated outgroup (which starts from the same consensus sequence and then considers subsequent fixations within the viral population). 

The steps for estimating the rate of adaptation within the genome of a virus are:
1. Setup and run a Nextstrain build
2. Compile a configuration file in `config/` that supplies the necessary locations of input data as well as additional biological metadata about the virus (such as it's cellular receptor or mode of transmission). The file `config/minimal_example.json` shows all of the necessary keys that must be present in the config file in order to calculate rates of adaptation.
3. Run the [`rate_of_adaptation.ipynb`](https://github.com/blab/adaptive-evolution/blob/master/adaptive-evolution-analysis/rate_of_adaptation.ipynb) notebook to compute the rate with an updating outgroup, or [`rate_of_adaptation_bhatt.ipynb`](https://github.com/blab/adaptive-evolution/blob/master/adaptive-evolution-analysis/rate_of_adaptation_bhatt.ipynb) to compute the rate with a constant outgroup.

Figures in the manuscripts can be generated using the following notebooks:
1. [Figure 1](https://github.com/blab/adaptive-evolution/blob/master/adaptive-evolution-analysis/Figure1.ipynb)
2. [Figure 2](https://github.com/blab/adaptive-evolution/blob/master/adaptive-evolution-analysis/Figure2.ipynb)
3. [Figure 3](https://github.com/blab/adaptive-evolution/blob/master/adaptive-evolution-analysis/Figure3.ipynb)
4. [Figure 4](https://github.com/blab/adaptive-evolution/blob/master/adaptive-evolution-analysis/Figure4.ipynb)
5. [Figure 5](https://github.com/blab/adaptive-evolution/blob/master/adaptive-evolution-analysis/Figure5.ipynb)
