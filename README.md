# **RECOGNICER  v1.0**

A coarse-graining based method for identifying multi-scale broad peaks from ChIP-seq data

## 1. Introduction

For details of the algorithm, please see

"*RECOGNICER: A coarse-graining approach for identifying broad domains from ChIP- seq data.*" Chongzhi Zang, Yiren Wang, and Weiqun Peng.(2020) doi: available soon...

This package is provided under the BSD-2-Clause lisense.


## 2. Installation

RECOGNICER v1.0 is implemented in python2. 
1. Prerequisites include `scipy` package. 

2. To install RECOGNICER, download the source package and setup the configurations in the script `RECOGNICER.sh`: set the variable `SRC_DIR` to be the directory where the uncompressed package is stored, e.g.,
`SRC_DIR=/home/packages/recognicer`

3. Make sure environment variables so that python modules can be executed smoothly.

4. Current supported genomes include: human hg38, hg19, hg18; mouse mm10, mm9, mm8; rat rn4; drosophila dm2, dm3; yeast sacCer1, pombe; and Arabidopsis tair8. For adding more customized species or genome versions, please see Additional notes below.


## 3. Running RECOGNICER

## 4. Output interpretation


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