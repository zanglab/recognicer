# **RECOGNICER v1.0**

A coarse-graining based method for identifying multi-scale broad peaks from ChIP-seq data

## 1. Introduction

## 2. Installation

## 3. Running RECOGNICER

## 4. Output notes


[https://genome.ucsc.edu/FAQ/FAQformat.html#format13](https://genome.ucsc.edu/FAQ/FAQformat.html#format13)

## 5. Additional notes

### Adding species/genome

Information of supported genomes is stored in `GenomeData.py`. If the species or genome version is not listed there, you can add it by editing `GenomeData.py` as follows:
1. Add a list of chromosomes as `${SPECIES}_chrom`
2. Add a dictionary of the length of each chromosomes as `${SPECIES}_chrom_lengths`
3. Modify the dictionary `species_chrom` by appending an element directing the genome name to the chrom list.
4. Modify the dictionary `species_chrom_lengths` by appending an element directing the genome name to the `chrom_lengths` dictionary.
5. If you do not want chrM in your analysis, simply delete `chrM` from the entries in `GenomeData.py`.

## 6. Contact

Email: zang@virginia.edu

[zanglab.org](http://zanglab.org)