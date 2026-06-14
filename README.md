# MHfinder

**MHfinder** is an open-source Python/R workflow for genome-wide microhaplotype identification and ancestry-informative panel construction from population genomic data. It was developed to support reproducible microhaplotype marker-panel design from PLINK-formatted datasets such as the 1000 Genomes Project.

The pipeline integrates:

1. population-aware SNP filtering by minor allele frequency (MAF),
2. LD-constrained microhaplotype construction,
3. stepwise conditional genome-wide association screening (SCGWAS),
4. MH phasing and summary statistics,
5. ancestry-prediction evaluation, and
6. marker-loss simulation for reduced-panel robustness assessment.

> **Manuscript context.** The accompanying study uses MHfinder to screen 1000 Genomes Project data and derive a compact ancestry-informative 48-MH panel. The same commands can be adapted to other population labels, subpopulation analyses, or targeted sequencing panel-design tasks.

---

## Repository layout

```text
core/                         Core Python modules
  snp_filter.py               Population-aware MAF filtering
  prep_MHs.py                 LD-constrained candidate MH generation
  scgwas_screen.py            Stepwise conditional GWAS screening
  utils.py                    Shared I/O and command helpers
extra/                        Optional analysis and plotting tools
  *.py                        LOOCV, validation, marker-loss and table scripts
  rscript/*.r                 R plotting and phasing scripts
examples/toy/                 Small public toy files for format demonstration
docs/                         Method and reproducibility notes
tests/                        Lightweight tests that do not require PLINK/R
```

---

## Requirements

### Required for the core Python package

- Python >= 3.8
- pandas >= 2.0.0
- openpyxl >= 3.1.0

### Required to run the full genomic workflow

- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- R >= 4.3
- R package `haplo.stats`

```r
install.packages("haplo.stats")
```

### Required for optional evaluation and plotting scripts

Python:

- numpy
- scipy
- scikit-learn
- statsmodels

R:

- optparse
- ggplot2
- readxl
- purrr
- readr
- dplyr
- patchwork
- cvms
- caret
- ggnewscale
- ggimage
- ggrepel
- adegenet
- plot3D

For convenience, see `environment.yml` for a conda-style environment specification.

---

## Installation

Clone the repository and install the Python package:

```bash
git clone https://github.com/Fun-Gene/MHfinder.git
cd MHfinder
pip install -e .
```

For optional analysis scripts:

```bash
pip install -e .[analysis]
```

Installation no longer requires PLINK or R to be present. Those executables are checked when commands that need them are run.

---

## Quick command-line workflow

### 1. Filter SNPs by population MAF

```bash
snp-filter-based-pop \
  --input-prefix data/1kg_phase3 \
  --sample-info sample_info.tsv \
  --maf-threshold 0.1 \
  --filter-mode union \
  --output-prefix results/filtered_snps
```

Filtering modes:

- `all`: MAF threshold across all samples
- `intersection`: MAF >= threshold in every population
- `union`: MAF >= threshold in at least one population

### 2. Generate candidate microhaplotypes

```bash
prep-MHs \
  --input-prefix data/1kg_phase3 \
  --snplist results/filtered_snps.snplist \
  --length-threshold 150 \
  --min-site-number 3 \
  --r2-threshold 0.8 \
  --output results/MHs_len150_n3_r0.8.txt
```

The output contains candidate MHs in the format:

```text
chr|SNP1-SNP2-...-SNPn
```

### 3. Screen MHs by stepwise conditional GWAS

```bash
screen-MHs-by-SCGWAS \
  --input-prefix data/1kg_phase3 \
  --sample-info sample_info.tsv \
  --mhs-file results/MHs_len150_n3_r0.8.txt \
  --p-threshold 5e-8 \
  --Ae-threshold 3.0 \
  --output results/Selected_MHs_5e-8_Ae3.txt
```

The selected-MH output uses three lines per candidate:

```text
>population_iteration_candidate
*indexSNP chromosome Ae(population)=value ...
SNP1 SNP2 SNP3 ...
```

---

## Optional downstream analysis tools

The `extra/` folder provides scripts used to generate evaluation tables and figures in the manuscript.

| Command | Purpose |
|---|---|
| `MLE_for_MHs_genotypes` | Estimate phased MH haplotypes using `haplo.stats` |
| `Prep_file_for_Infocalc` | Convert MH genotypes to Infocalc input |
| `Get_MHlist_with_max_In` | Select one MH with maximum Rosenberg's informativeness per duplicated index-SNP group |
| `Filter_MH_genotype_file_based_MHlist` | Filter an MH genotype table by selected MH list |
| `LOOCV_of_MHs` | Leave-one-out Naive Bayes ancestry-prediction evaluation |
| `Naive_Bayes_random_kept_MHs` | Simulate marker loss by random retained-marker subsets |
| `Validation_of_MHs_using_test_set` | Evaluate trained MH model on an independent test set |
| `Draw_manhattan_qq_of_SCGWAS` | Draw Manhattan/Q-Q plots from SCGWAS results |
| `Draw_PCA_of_MHs` | Draw 3D PCA plots from MH genotypes |
| `AUC_ACC_decline_with_n_MHs` | Summarize AUC/accuracy decline under marker loss |
| `Stat_MHs_indexSNPs_info` | Summarize selected index-SNP frequencies and regression statistics |
| `Organize_MHs_information_table` | Organize selected MH annotations and forensic parameters |
| `Stat_allele_freq_of_MHs` | Calculate MH allele-frequency tables |
| `Calculate_Ae_of_MHs` | Calculate effective number of alleles across populations |

The legacy misspelled command `Organize_MHs_infomation_table` is retained as an alias for backward compatibility.

---

## Input file formats

### Sample information table

A tabular file in `.xlsx`, `.csv`, `.tsv`, or `.txt` format. The first three columns are interpreted as:

```text
FamilyID    SampleID    Population
```

Example:

```text
FamilyID    SampleID    Population
FAM001      HG00001     EUR
FAM002      HG00002     AFR
```

### PLINK files

The core workflow expects a standard PLINK binary prefix:

```text
<prefix>.bed
<prefix>.bim
<prefix>.fam
```

---

## Examples

The `examples/toy/` directory contains small public demonstration files for input/output format inspection and parser tests. These files are synthetic and are not intended to reproduce the manuscript's 1000 Genomes results.

A full 1000 Genomes reproduction guide is provided in `docs/reproduce_1000genomes_workflow.md`.

---

## Testing

Lightweight tests that do not require PLINK or R can be run with:

```bash
python -m unittest discover -s tests
```

These tests check importability, table parsing, CLI help messages, and candidate-MH parsing behavior. Full end-to-end genomic execution requires PLINK/R and a PLINK dataset.

---

## Reproducibility notes

The manuscript analysis used public 1000 Genomes Project data. Because full WGS files are large, this repository provides:

- command-level workflow documentation,
- toy data for format demonstration,
- scripts used for evaluation and plotting, and
- table outputs from the selected panel where appropriate.

For manuscript submission, we recommend archiving a tagged GitHub release on Zenodo and citing the resulting DOI.

---

## Citation

If you use MHfinder, please cite the associated manuscript and the archived software release. A machine-readable citation template is provided in `CITATION.cff`.

---

## License

MHfinder is released under the GNU General Public License v3.0. See `LICENSE`.
