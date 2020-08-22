# **RECOGNICER  v1.0**

A coarse-graining based method for identifying multi-scale broad peaks from ChIP-seq data

## Introduction

For details of the algorithm, please see

"*RECOGNICER: A coarse-graining approach for identifying broad domains from ChIP-seq data.*" Chongzhi Zang, Yiren Wang, and Weiqun Peng.(2020) doi: available soon...

This package is provided under the [BSD-2-Clause](https://opensource.org/licenses/BSD-2-Clause) license. Please cite the above paper in your publication if you use this algorithm or package for analyzing the data in your work.


## Installation

RECOGNICER v1.0 is implemented in python2. 
1. Prerequisites include the `scipy` package. 

2. To install RECOGNICER, download the source package and setup the configurations in the script `RECOGNICER.sh`: set the variable `SRC_DIR` to be the directory where the uncompressed package is stored, e.g.,
`SRC_DIR=/home/packages/recognicer`

3. Make sure environment variables so that python modules can be executed smoothly.

4. Current supported genomes include: human hg38, hg19, hg18; mouse mm10, mm9, mm8; rat rn4; drosophila dm2, dm3; yeast sacCer1, pombe; and Arabidopsis tair8. For adding more customized species or genome versions, please see Additional notes below.


## Running RECOGNICER

The input of RECOGNICER is two ChIP-seq mapped read data files, for the ChIP signal sample and the control sample, respectively in [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format. To run RECOGNICER, please first set the configuration parameters in `RECOGNICER.sh`, including a customized sample name `SAMPLE` (default: Sample), species/genome version `SPECIES` (default: hg18), and output directory `OUT_DIR`. Then use the command

`sh RECOGNICER.sh ChIP.bed control.bed 0.01`

The 3 input arguments are ChIP data file, control data file, and FDR threshold. The whole process may take several hours. 

Although not required, you can also tune the following parameters as you like:
- `WINDOW_SIZE`: Resolution of the RECOGNICER result, in bp. default: 200 (mononucleosome+linker).
- `FRAGMENT_SIZE`: ChIP DNA fragment size, used for determination of the shift size from the 5' start location of a sequence read towards the center of the DNA fragment represented by the read. default: 150 (mononucleosome).
- `GENOME_FRACTION`: Effective mappable fraction of the whole genome, dependent on the read length. default: 0.74 for the human genome.
- `CHIPTHRESHOLD`: Number of copies of identical reads allowed in a dataset. default: 1.

It is strongly NOT recommended to change other parameters.

**Caution**: please do not run multiple instances of RECOGNICER scripts originated from the same directory in parallel. Each instance of RECOGNICER script generates temporary files with hard-coded names. If multiple instances of RECOGNICER scripts originated from the same directory are run in parallel, their temporary files would interfere with each other, resulting in absurd outcomes. If you need to run multiple instances of RECOGNICER simultaneously, please make sure to run them under separate directories to avoid interference.


## Output interpretation

The output of RECOGNICER includes the following data files:
- `$SAMPLE-nonredundant.bed`: Non-redundant reads filtered from the input ChIP dataset.
- `$SAMPLE-W200.bedgraph`: Sample read pile-up track file, window size (200bp) resolution, in [bedGraph](https://genome.ucsc.edu/FAQ/FAQformat.html#format1.8) format, for genome browser visualization.
- `$SAMPLE.cgsummary`: Summary of all candidate domains called from coarse-graining, with their complete information and statistical assessment, in the following format: chrom, start, end, ChIP read count, control read count, P-value, fold change, Q-value.
- `$SAMPLE-fdr0.01_broadPeak.bed`: Final identified significant domains, in [ENCODE broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13) format. The score in the 5th column is the fold change comparing ChIP and control read counts. The signalValue in the 7th column is the ChIP read count.


## Additional notes

### Adding customized species/genomes

Information of supported genomes is stored in `GenomeData.py`. If the species or genome version of your data is not included, you can add it by editing `GenomeData.py` as follows:
1. Add a list of chromosomes as `${SPECIES}_chrom`
2. Add a dictionary of the length of each chromosomes as `${SPECIES}_chrom_lengths`
3. Modify the dictionary `species_chrom` by appending an element directing the genome name to the chrom list.
4. Modify the dictionary `species_chrom_lengths` by appending an element directing the genome name to the `chrom_lengths` dictionary.
5. If you do not want chrM in your analysis, simply delete `chrM` from the entries in `GenomeData.py`.


## Contact

Email: zang@virginia.edu

[zanglab.org](http://zanglab.org)
