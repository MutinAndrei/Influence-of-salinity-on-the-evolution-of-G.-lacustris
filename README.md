# Halotolerant-gammarus-lacustris-phylogenetics

[![License: Custom](https://img.shields.io/badge/License-Custom-blue)](LICENSE.md)
[![Status](https://img.shields.io/badge/Status-Article%20Submitted-lightgrey)]()

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
| **This study** | Phylogenetic tree images in PNG and SVG formats used in this study | [Figure](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/results/Figure) |
| **This study** | Phylogenetic tree generation files used in this study | [trees](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/results/trees) |
| **This study** | BASH script for generating the phylogenetic tree, R script for retrieving FASTA sequence files from GenBank, R script for visualizing phylogenetic trees | [scripts](https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris/tree/main/scripts) |

## ⚙️ Reproducing the analysis

To reproduce the phylogenetic analysis presented in the article, the following steps must be performed.

### 1. **Cloning the repository**

```
git clone https://github.com/MutinAndrei/Influence-of-salinity-on-the-evolution-of-G.-lacustris.git
cd Influence-of-salinity-on-the-evolution-of-G.-lacustris
```

### 2. **Installing dependencies**

The following software was used for the phylogenetic analysis:

- IQ-TREE (v2.0.7) — The program can be installed via the link [IQ-TREE](https://iqtree.github.io/)

- MAFFT (v7.511) — The program can be installed via the link [MAFFT](https://mafft.cbrc.jp/alignment/software/linux.html)
 
- TrimAl (v1.5) — The program can be installed via the link [TrimAl](https://github.com/inab/trimal)

For tree visualization, the R programming language (v4.5.0) and the following packages were used:

- tidyverse (v2.0.0)
- ggtree (v3.16.0)
- ape (v5.8-1)
- ggplot2 (v4.0.1)
- dplyr (v1.1.4)
- stringr (v1.5.1)


### 3. **Running commands**

3.1. **Cycle MAFFT alignment**

Navigate to the `/data/raw_seq` folder and run the following command:

```
for x in *.fa; do mafft $x > $x.aln; done
```

The execution of this command will result in the creation of `.aln` files in the `raw_seq` folder.

3.2. **Trim alignment cycle**

After that, trim the alignment using the following command (specifying the path to the TrimAl program beforehand):

```
for x in *.aln; do /path/to/trimal -in $x -out $x.trim -automated1; done
```

The execution of this command will result in the creation of files containing trimmed sequences with the `.trim` extension.

3.3. **Creating a folder with trimmed files** 

Create a directory named Trim_file

Move the resulting `.trim` files into the directory named Trim_file
	
Next, navigate to the `partitions` folder and copy the file `example.nex` into the `raw_seq` directory.

3.4. **Building a tree**

To build the phylogenetic tree, run the command:

```
iqtree2 -p Trim_file/ -p 8genes_partitions.nex -B 1000
```

The main result of executing this command will be the generation of a number of tree files, the primary one being the file with the `.treefile` extension.

### 4. Tree visualization

4.1. To visualize the trees, copy the newly generated file with the `.treefile` extension into the `scripts` folder.

4.2. Run the R script named R_script_for_rendering_trees.R

4.3. The result of executing the script will be phylogenetic tree images in PNG and SVG formats.

### 5. Statistical analysis

5.1. To perform the statistical analysis, take the scripts `IBD_script.R` and `Blomberg_K_and_Pagel.R` from the `scripts` folder and place them in the same directory as the `.treefile` tree file.

## 📬 Contacts

For any questions regarding the data or analysis, please contact:

Ekaterina Shchapova - shchapova.katerina@gmail.com

Andrei Mutin - andreimutin97@gmail.com
