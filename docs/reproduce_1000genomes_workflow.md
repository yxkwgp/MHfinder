# Reproducing the 1000 Genomes MHfinder workflow

This document records the command-level workflow used for a full-scale MHfinder analysis on public 1000 Genomes Project data. It is intended as a reproducibility guide for the manuscript, not as a small tutorial; full WGS processing requires substantial disk space and compute time.

## 1. Data sources

- 1000 Genomes Project Phase 3 whole-genome variants for panel development.
- Expanded high-coverage 1000 Genomes release for independent validation.
- Population/sample metadata from the 1000 Genomes Project.

Record exact URLs, release dates, checksums, and genome build in the final manuscript submission package.

## 2. Convert VCF to PLINK

Example command pattern per chromosome:

```bash
vcftools \
  --vcf ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  --plink \
  --out plink/chr${CHR}
```

Merge chromosome-level PLINK files into a genome-wide binary prefix following PLINK documentation. Exclude sites with both SNP and indel records before merging.

## 3. Prepare sample information

Create a three-column file:

```text
FamilyID    SampleID    Population
```

For the manuscript analysis, `Population` is one of:

```text
AFR, AMR, EAS, EUR, SAS
```

## 4. Run MHfinder core workflow

```bash
mkdir -p results

snp-filter-based-pop \
  --input-prefix plink/1kg_phase3_merged \
  --sample-info metadata/sample_info.tsv \
  --maf-threshold 0.1 \
  --filter-mode union \
  --output-prefix results/filtered_snps

prep-MHs \
  --input-prefix plink/1kg_phase3_merged \
  --snplist results/filtered_snps.snplist \
  --length-threshold 150 \
  --min-site-number 3 \
  --r2-threshold 0.8 \
  --output results/MHs_len150_n3_r0.8.txt

screen-MHs-by-SCGWAS \
  --input-prefix plink/1kg_phase3_merged \
  --sample-info metadata/sample_info.tsv \
  --mhs-file results/MHs_len150_n3_r0.8.txt \
  --p-threshold 5e-8 \
  --Ae-threshold 3.0 \
  --output results/Selected_MHs_5e-8_Ae3.txt
```

## 5. Downstream evaluation

After SCGWAS selection, use the `extra/` tools to phase MHs, compute allele frequencies/Ae, run LOOCV, simulate marker loss, and validate on the independent test set. Example command sequence:

```bash
MLE_for_MHs_genotypes \
  -p plink/1kg_phase3_merged \
  -m results/Selected_MHs_5e-8_Ae3.txt \
  -o results/MHs_genotypes.tsv

LOOCV_of_MHs \
  -g results/MHs_genotypes.tsv \
  -s metadata/sample_info.tsv \
  -o results/loocv_48MHs

Naive_Bayes_random_kept_MHs \
  -g results/MHs_genotypes.tsv \
  -s metadata/sample_info.tsv \
  -r 2026 \
  -t 10000 \
  -d results/random_marker_loss
```

## 6. Reporting checklist for manuscript reproducibility

Before submission, record:

- data URLs and checksums,
- software versions for PLINK, R, haplo.stats, Python and MHfinder,
- exact command lines,
- CPU/RAM used,
- wall-clock runtime for each stage,
- random seeds,
- output hashes for key result files.
