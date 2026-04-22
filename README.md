# Halotolerant-gammarus-lacustris-phylogenetics

**Repository for the article:**
> *"Halotolerant Gammarus lacustris Sars, 1863 (Amphipoda) in water bodies of the Baikal region with different salinity regimes: role in ecosystems and perspectives in aquaculture"* 

> The title and link to the article will be here (The article is currently in the publication process; this information will be updated once the journal issue with the article is released).

## 📌 Overview
This repository contains data and scripts for the phylogenetic analysis described in the article "Halotolerant Gammarus lacustris Sars, 1863 (Amphipoda) in water bodies of the Baikal region with different salinity regimes: role in ecosystems and perspectives in aquaculture". All data stored in this repository are necessary for the complete reproduction of the multilocus phylogenetic analysis presented in the manuscript.

To reconstruct the evolutionary relationships among populations of *G. lacustris* inhabiting water bodies with contrasting salinity, a dataset of eight molecular markers following [**Hou et al. (2022)**](https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16160) and [**Gurkov et. al. (2019)**](https://doi.org/10.1186/s12862-019-1470-8):

- Mitochondrial genes: *COI*, *COIII*, *12S*, *16S*;
- Nuclear genes: *AK*, *EF1α*, *HSP70*, *PEPCK*, *28S*.

In addition to the sequences from the studies by [**Hou et al. (2022)**](https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16160) and [**Gurkov et al. (2019)**](https://doi.org/10.1186/s12862-019-1470-8) we obtained our own *COI* sequences of *G. lacustris* from water bodies of the Baikal region and included them in the phylogenetic analysis.

The repository includes:
- **Original sequences:** FASTA files for each gene prior to the alignment step.
- **Sample metadata:** A table containing information on location, salinity, and source for each sequence.
-- **Partition and model file:** A Nexus file that defines the partitioning scheme for eight genes and the **nucleotide substitution models from Hou et al. (2022)**.
- **Phylogenetic trees:** The final tree in Newick format and the complete IQ-TREE log file.
- **Scripts:** R code for data preparation and visualization, as well as a bash script for running the entire computational pipeline.

## 🧬 Data Availability

| Resource | Description | Source |
| :--- | :--- | :--- |
| **This study** | Original chromatograms (`.ab1`) | Available upon request |
| **This study** | Consensus sequences | NCBI GenBank: `[Accession numbers pending]` |
| **Hou et al. (2022)** | Sequences of 8 genes of the *Gammarus lacustris* clade  | [Hou et al. (2022)](https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16160) |
| **Gurkov et al. (2019)** | COI sequences of *Gammarus lacustris* from GenBank | [Gurkov et al. (2019)](https://doi.org/10.1186/s12862-019-1470-8) |
| **This study** | Table with data on sampling locations and water body salinity | [metadata](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/data/metadata) |
| **This study** | Partition file with substitution models from [**Hou et al. (2022)**](https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16160) | [partitions](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/data/partitions) |
| **This study** | FASTA files for each gene containing sequences from [**Hou et al. (2022)**](https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.16160), [**Gurkov et al. (2019)**](https://doi.org/10.1186/s12862-019-1470-8), and this study | [raw_seq](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/data/raw_seq) |
| **This study** | Final aligned and trimmed sequences | [alignments](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/results/alignments) |
| **This study** | Phylogenetic tree images in PNG and SVG formats used in this study | [figures](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/results/figures) |
| **This study** | Phylogenetic tree generation files used in this study | [trees](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/results/trees) |
| **This study** | BASH script for generating the phylogenetic tree, R script for retrieving FASTA sequence files from GenBank, R script for visualizing phylogenetic trees | [scripts](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/scripts) |

## 📬 Contacts

For any questions regarding the data or analysis, please contact:

Shchapova Ekaterina Pavlovna - shchapova.katerina@gmail.com

Mutin Andrei - andreimutin97@gmail.com
