# MHfinder method summary

MHfinder implements a population-aware computational workflow for discovering ancestry-informative microhaplotypes (MHs) from genome-wide genotype data.

## Stage 1: population-aware SNP filtering

Input SNPs are filtered by minor allele frequency (MAF). In addition to global filtering, MHfinder supports population-aware union and intersection filters so that variants informative in one population are not discarded solely because they are rare globally.

## Stage 2: LD-constrained MH generation

Candidate MHs are generated from SNPs located within a user-defined physical window, for example 150 bp. A candidate must contain at least a user-defined number of SNPs and must satisfy an LD constraint. The implementation uses PLINK pairwise LD output and a Bron-Kerbosch clique search to identify maximal SNP combinations meeting the pairwise LD criterion.

## Stage 3: stepwise conditional GWAS screening

For each population label, MHfinder encodes membership as a one-vs-rest phenotype and runs PLINK logistic regression. At each iteration, the most significant SNP that maps to at least one eligible MH is selected. Previously selected SNPs are added as covariates in the next iteration, and SNPs from selected MHs are excluded from subsequent rounds. The process stops when no remaining SNP exceeds the p-value threshold or an optional maximum number of iterations is reached.

## Stage 4: optional downstream evaluation

The optional scripts in `extra/` support:

- MH haplotype phasing,
- Rosenberg's informativeness input preparation,
- effective allele number calculation,
- leave-one-out cross-validation,
- independent test-set validation,
- random marker-loss simulation,
- PCA/Manhattan/Q-Q/confusion-matrix plotting.

## Recommended manuscript reporting

For a software or methods manuscript, report:

- input data release and genome build,
- sample labels and sample sizes,
- MAF, length, SNP-count, LD, p-value and Ae thresholds,
- number of SNPs and candidate MHs after each stage,
- runtime and peak memory for each stage,
- comparison against alternative marker-ranking strategies,
- random seeds and software versions.
