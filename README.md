# Tumor Phylogeny Simulator

Create a phylogenetic simulator to generate artificial data sets generated using a variety of parameters influencing the complexity and progression of a tumor.

* [Slides](slides.pdf)
* [Report](report.pdf)

## Problem

Tumors are naturally heterogeneous, meaning each individual tumor is comprised of many cell types each with theirown genetic and epigenetic profile.  When patients undergo cancer treatment, tumor samples can be taken and sub-sequently sequenced.  The most common sequencing technique is called bulk sequencing where the sequenced readscan align to any of the cell types from the original tumor.  A majority of clinics use bulk sequencing over single cellsequencing as single cell sequencing is more expensive and time consuming.

As a consequence, there exist many tools designed to \unmix" or \deconvolve" these bulk sequences into their original cell types.  [6][5][1][2][9][3][4] Some of these approaches use single nucleotide polymorphisms (SNPs) [4], some use copy number variations (CNVs) [9], while others use both [6][5][1][3].  Each of these tools take their own custom input and output as well as simulating their own data making it difficult to adequately compare tools.  Previous attempts to create generalized simulators have not covered all possible mutation events and can therefore not be used to generate a standardized test data set.  One simulator, PyVolve, uses a continuous time Markov model to build a tumor phylogeny with SNP mutations but fails to infer CNVs [8].  Another called Pomegranate infers SNPs and CNVs but does not maintain a genome so each CNV acts independently [7].  There is also no information produced on structural variants and methylation.

## Solution

To  create  a  standardized  data  set  for  the  deconvolution  tools,  we  propose  creating  a  phylogenetic  simulator  that will use a continuous time markov model to determine when and which mutations occur next at each time.  When mutations occur, a branch in the phylogeny will be created bearing the same genetic profile as it's parent but now containing  this  additional  mutation.   The  set  of  mutations  we  aim  to  create  will  be:
1. CNV inversion 
2. CNV amplification
3. CNV deletion

Each node in the tree will keep track of the exact genetic profile of a cell type and the total numberof cells with that genetic profile.  We will then artificially "mix" these cells to produce mixed variant allele frequencies(VAFs) for SNPs, mixed copy number values for CNVs, mixed breakpoint copy numbers for SVs, and mixed variant allele frequencies for methylated CpGs.  Our final task will be to use this simulator to create a standardized data setwith many simulated samples binned into simulated "patients" to be used by the tools above.  We will adhere to thevariant calling format (VCF) as is standard for documenting variants.

## What did we do:
We simulated 288 data sets generated using a variety of parameters influencing the complexity and progression of a tumor. All data can be found at https://cmu.box.com/s/thq5ba255f24av575ax5xhk9kvp2syve. The parameters that were used for constructing this data were: number of samples, number of mutations, mean mutation length, value choices for parameters for beta distributions, and choice in variance for generating a mixed sample. Along with these parameters, we also created a small data set labeled "toy" and larger data set labeled "real". The toy data set contains five chromosomes of lengths $10000, 15000, 20000, 25000, 30000$ while the real data set contains the 22 human chromosomes excluding X and Y chromosomes.


## Tasks

- [x] make github
- [x] make GeneProf
- [x] CNV duplication
- [x] CNV inversion
- [x] CNV deletion
- [x] Research and understand CTMMs - all
- [x] Find the probability/rate of future mutations based on past mutations - shef
- [x] Tree generation where each node is a GeneProf - murtaza, jesse
- [x] Tree Collapse (combine nodes) and make C - jesse, shef
- [x] Generate U and mix U * C --> F - jesse, murtaza
- [ ] VCF formatize - future work
- [x] Write report - all 
- [x] Make ppt - all
