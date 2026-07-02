# MHfinder - a stepwise conditional GWAS pipeline for genome-wide discovery of ancestry-informative microhaplotypes  
  
**MHfinder** is an Python module for genome-wide microhaplotype identification and ancestry-informative panel construction from population genomic data.  
  
The module integrates:  
1. population-aware SNP filtering by minor allele frequency (MAF)  
2. genome-wide microhaplotype construction  
3. candidate ancestry informative microhaplotypes finding by stepwise conditional genome-wide association screening (SCGWAS)  
4. optional analysis script for evaluation of MH-panel  
  
## Requirements  
  
### Required to run the full genomic workflow  
  
- `Python` >= 3.8  
- Python module `pandas` >= 2.0.0  
- Python module `openpyxl` >= 3.1.0  
- `R` >= 4.3  
- R package `haplo.stats`  
- PLINK 1.9 (https://www.cog-genomics.org/plink/)  
  
### Required for optional analysis scripts  
  
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
  
## Installation  
  
### 1. Install PLINK 1.9  
Download PLINK 1.9 (https://www.cog-genomics.org/plink/)  
Add to PATH  
Verify installation  
```bash  
plink --version  
# Should display: PLINK v1.9...  
```  
  
### 2. Install R (>= 4.3)  
Windows: Download from https://cran.r-project.org/bin/windows/base/  
macOS: Download from https://cran.r-project.org/bin/macosx/  
Linux: Download from https://cran.r-project.org/bin/linux/  
```bash  
# Verify installation  
Rscript --version  
# Should display: R version 4.3.x or higher  
```  
  
### 3. Install Required R Package  
```bash  
# Start R  
install.packages("haplo.stats")  
```  
  
### 4. Download and Install Python Module  
Clone the repository  
```bash  
git clone https://github.com/Fun-Gene/MHfinder.git  
```  
  
Install the Python module without extra analysis scripts:  
```bash  
cd MHfinder  
pip install -e .  
```  
  
Install the Python module with extra analysis scripts:  
```bash  
cd MHfinder  
pip install -e .[analysis]  
```  
  
## Quick command-line workflow  
The input files for this workflow are the PLINK binary files (.bed, .bim, .fam) and a sample information file that records the sample metadata for this set of PLINK files.  
The sample info file should be a tabular file in `.xlsx`, `.csv`, `.tsv`, or `.txt` format.  
Example:  
```text  
FID IID population  
FAM001  HG00001 EUR  
FAM002  HG00002 AFR  
```  
Note: The FID and IID columns in the sample information file must match the corresponding columns in the PLINK .fam file.  
  
### 1. Filter SNPs by population MAF  
  
```bash  
snp-filter-based-pop \  
  --input-prefix data/1kg_phase3 \  
  --sample-info sample_info.tsv \  
  --maf-threshold 0.1 \  
  --filter-mode union \  
  --output-prefix results/filtered_snps  
```  
  
**Parameter Description:**  
- `input_prefix`: PLINK binary file prefix  
- `sample_info`: Sample information file  
- `maf_threshold`: MAF threshold  
- `filter_mode`: Filtering mode  
  - `'all'`: MAF across all samples  
  - `'intersection'`: MAF ≥ threshold in all populations  
  - `'union'`: MAF ≥ threshold in at least one population  
- `output_prefix`: Output file prefix  
  
The output is a snplist file contains SNPs passing filter  
  
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
**Parameter Description:**  
- `input_prefix`: PLINK binary file prefix  
- `snplist`: Filtered snplist file  
- `length_threshold`: Length threshold (distance threshold between two endpoint SNPs of microhaplotype)  
- `min_site_number`: Minimum number of SNPs in microhaplotype  
- `r2_threshold`: Strong linkage threshold (r2 between any two SNPs in microhaplotype must be below this value)  
- `output`: Final output filename  
  
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
  
**Parameter Description:**  
- `input_prefix`: PLINK binary file prefix  
- `sample_info`: Sample information file  
- `mhs_file`: Microhaplotype file from previous step  
- `p_threshold`: GWAS p-value threshold  
- `Ae_threshold`: minimum average effective number of alleles (Ae) threshold  
- `loop_times`: Maximum GWAS loop times (optional, default unlimited)  
- `output`: Output file name  
  
The selected-MH output uses three lines per candidate:  
  
```text  
>population_iteration_candidate  
*indexSNP chromosome Ae(population)=value ...  
SNP1 SNP2 SNP3 ...  
```  
  
## Optional downstream analysis tools  
  
If this module was installed in [analysis] mode, the following command is available to produce a series of assessment tables and plots for the AIMs that were filtered above. Add `--help` to the command for usage information.  
  
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
  
## Demo data  
  
The demo.zip contains the input files we prepared for trying out this Python module, including a set of PLINK binary files (subset from chromosome 22 of the 1000 Genomes Phase 3 data, retaining only SNPs) and a matched sample information table. You can try applying the module—from start to obtaining candidate markers, it typically takes only about 10 minutes on a PC.  
  
## Testing  
  
Lightweight tests that do not require PLINK or R can be run with:  
  
```bash  
python -m unittest discover -s tests  
```  
These tests check importability, table parsing, CLI help messages, and candidate-MH parsing behavior.  
  
## Notes  
  
1. Ensure PLINK 1.9 is properly installed and accessible from command line  
2. Ensure R 4.3+ is properly installed and accessible from command line  
3. Ensure R package 'haplo.stats' is installed  
4. Input file paths must exist and have correct format  
5. Temporary files will be generated in current directory, recommend running in dedicated working directory  
