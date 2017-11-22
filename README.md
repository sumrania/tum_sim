# Tumor Phylogeny Simulator

## Problem

Tumors are naturally heterogeneous, meaning each individual tumor is comprised of many cell types each with theirown genetic and epigenetic profile.  When patients undergo cancer treatment, tumor samples can be taken and sub-sequently sequenced.  The most common sequencing technique is called bulk sequencing where the sequenced readscan align to any of the cell types from the original tumor.  A majority of clinics use bulk sequencing over single cellsequencing as single cell sequencing is more expensive and time consuming.

As a consequence, there exist many tools designed to \unmix" or \deconvolve" these bulk sequences into theiroriginal cell types.  [6][5][1][2][9][3][4] Some of these approaches use single nucleotide polymorphisms (SNPs) [4], someuse copy number variations (CNVs) [9], while others use both [6][5][1][3].  Each of these tools take their own custominput and output as well as simulating their own data making it difficult to adequately compare tools.  Previous attempts to create generalized simulators have not covered all possible mutation events and can therefore not be usedto generate a standardized test data set.  One simulator, PyVolve, uses a continuous time Markov model to builda tumor phylogeny with SNP mutations but fails to infer CNVs [8].  Another called Pomegranate infers SNPs andCNVs but does not maintain a genome so each CNV acts independently [7].  There is also no information producedon structural variants and methylation.

## Solution

To  create  a  standardized  data  set  for  the  deconvolution  tools,  we  propose  creating  a  phylogenetic  simulator  thatwill use a continuous time markov model to determine when and which mutations occur next at each time.  Whenmutations occur, a branch in the phylogeny will be created bearing the same genetic profile as it's parent but nowcontaining  this  additional  mutation.   The  set  of  mutations  we  aim  to  create  will  be:
1. SNP  transition
2. SNP transversion
3. CNV amplification
4. CNV deletion
5. Structural variant (SV) inversion
6. SV translocation
7. CpG methylation.  

Each node in the tree will keep track of the exact genetic prole of a cell type and the total numberof cells with that genetic profile.  We will then articially "mix" these cells to produce mixed variant allele frequencies(VAFs) for SNPs, mixed copy number values for CNVs, mixed breakpoint copy numbers for SVs, and mixed variantallele frequencies for methylated CpGs.  Our final task will be to use this simulator to create a standardized data setwith many simulated samples binned into simulated "patients" to be used by the tools above.  We will adhere to thevariant calling format (VCF) as is standard for documenting variants.
