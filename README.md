# MHfinder - Finding candidate ancestry inference microhaplotypes

## Project Overview

This is a Python module for finding candidate ancestry inference microhaplotypes, integrating three core functions:
1. **Population-based MAF SNP filtering** - Filter SNPs based on Minor Allele Frequency (MAF)
2. **Microhaplotype generation** - Generate microhaplotypes based on filtered SNPs
3. **Stepwise conditional GWAS-based microhaplotype screening** - Screen candidate microhaplotypes using stepwise conditional GWAS method

## Requirements

- **PLINK 1.9**: PLINK 1.9 version must be installed, check version by running `plink --version`
- **R >= 4.3**: R version 4.3 or higher must be installed
- **Python >= 3.8**: Python version 3.8 or higher required

### R Package Requirements
- **haplo.stats**
  ```r
  install.packages("haplo.stats")
  ```

### Python Package Requirements
- **pandas >= 2.0.0**
- **openpyxl >= 3.1.0**

## Installation Instructions

### 1. Install PLINK 1.9
```bash
# Download PLINK 1.9
# Visit https://www.cog-genomics.org/plink/
# Download version suitable for your system

# Extract and add to PATH
# Windows: Add plink.exe directory to system PATH
# Linux/Mac: Copy plink executable to /usr/local/bin/ or add to PATH

# Verify installation
plink --version
# Should display: PLINK v1.9...
```

### 2. Install R (>= 4.3)
```bash
# Windows: Download from https://cran.r-project.org/bin/windows/base/
# macOS: Download from https://cran.r-project.org/bin/macosx/
# Linux: Download from https://cran.r-project.org/bin/linux/

# Verify installation
Rscript --version
# Should display: R version 4.3.x or higher
```

### 3. Install Required R Package
```r
# Start R and install the required package
R
install.packages("haplo.stats")
```

### 4. Download and Install Python Module
```bash
# Download the repository from GitHub
# 1. Click the green "Code" button
# 2. Select "Download ZIP"
# 3. Extract the downloaded ZIP file to your desired location
# 4. Navigate to the extracted folder
cd MHfinder
# Install module (will automatically check all dependencies)
pip install -e .
# If any dependency version doesn't meet requirements, installation will fail and display error message
```

## Module Functions

### 1. Population-based MAF SNP Filtering (filter_snps_by_maf)

Filter SNPs based on Minor Allele Frequency (MAF), supporting multiple filtering modes.

**Parameter Description:**
- `input_prefix`: PLINK binary file prefix
- `sample_info`: Sample information file
- `maf_threshold`: MAF threshold
- `filter_mode`: Filtering mode
  - `'all'`: MAF across all samples
  - `'intersection'`: MAF ≥ threshold in all populations
  - `'union'`: MAF ≥ threshold in at least one population
- `output_prefix`: Output file prefix

**Returns:** Generated snplist file path

### 2. Microhaplotype preparation (prep-MHs)

Prepare microhaplotypes based on filtered SNPs.

**Parameter Description:**
- `input_prefix`: PLINK binary file prefix
- `snplist`: Filtered snplist file
- `length_threshold`: Length threshold (distance threshold between two endpoint SNPs of microhaplotype)
- `min_site_number`: Minimum number of SNPs in microhaplotype
- `r2_threshold`: Strong linkage threshold (r2 between any two SNPs in microhaplotype must be below this value)
- `output`: Final output filename

**Returns:** Generated microhaplotype file path

### 3. Stepwise Conditional GWAS-based Microhaplotype Screening (screen_mhs_by_scgwas)

Screen candidate microhaplotypes using stepwise conditional GWAS method.

**Parameter Description:**
- `input_prefix`: PLINK binary file prefix
- `sample_info`: Sample information file
- `mhs_file`: Microhaplotype file from previous step
- `p_threshold`: GWAS p-value threshold
- `Ae_threshold`: minimum average effective number of alleles (Ae) threshold
- `loop_times`: Maximum GWAS loop times (optional, default unlimited)
- `output`: Output file name

**Returns:** Screened microhaplotype file path

## Usage

### Command Line Usage

After installation, you can use three independent commands:

```bash
# 1. Filter SNPs based on population MAF
snp-filter-based-pop \
    --input-prefix / -i <plink_prefix> \
    --sample-info / -s sample_info.csv \
    --maf-threshold / -t 0.1 \
    --filter-mode / -m union \
    --output-prefix / -o filtered_snps

# 2. Generate microhaplotypes
prep-MHs \
    --input-prefix / -i <plink_prefix> \
    --snplist / -s filtered_snps.snplist \
    --length-threshold / -l 150 \
    --min-site-number / -n 3 \
    --r2-threshold / -r 0.8 \
    --output / -o MHs_len150_n3_r0.8.txt

# 3. Screen microhaplotypes based on stepwise conditional GWAS
screen-MHs-by-SCGWAS \
    --input-prefix / -i <plink_prefix> \
    --sample-info / -s sample_info.csv \
    --mhs-file / -m MHs_len150_n3_r0.8.txt \
    --p-threshold / -p 5e-8 \
    --Ae-threshold / -a 3.0 \
    --output / -o Selected_MHs_5e_8_Ae3.txt
```

## File Format Description

### Sample Information File Format
- Must contain columns: FamilyID, SampleID, Population
- Supported formats: xlsx, csv, tsv, txt

### Microhaplotype File Format
- Each line format: `chr|SNP1-SNP2...SNPn`
- Example: `1|rs123-rs456-rs789`

### Selected Ancestry Inference Microhaplotype File Format
- Each microhaplotype is represented by 3 lines:
  - **Line 1**: `>genotype_loop_index` - Contains the genotype and loop iteration information
  - **Line 2**: `*indexSNP\tchromosome\tDiversity=value\tAe(population)=value` - Contains key SNP, chromosome, diversity value, and effective number of alleles information
  - **Line 3**: `SNP1\tSNP2\tSNP3\t...` - Contains all SNPs in the microhaplotype, separated by tabs

**Example:**
```
>AFR_1_1
*rs123456	1	Diversity=0.75	Ae(AFR)=3.2	Ae(AMR)=2.8
rs123	rs456	rs789
```

## Extra Tools

If you would like to obtain the figures and tables presented in our paper, which include the evaluation, statistics, and visualization of the markers screened by our tool, we recommend that you run the following tools in order. They are optional scripts from the `extra/` folder.

### Dependencies for `extra` Tools

**Python**
- **sklearn**
- **scipy**
- **statsmodels**

**R**
- **optparse**
- **ggplot2**
- **readxl**
- **purrr**
- **readr**
- **dplyr**
- **patchwork**
- **cvms**
- **caret**
- **ggnewscale**
- **ggimage**
- **ggrepel**
- **adegenet**
- **plot3D**

### Tool Usage

#### `MLE_for_MHs_genotypes`
Estimate MH genotypes by maximum likelihood with R phasing support. The <mhs_info_file> is selected ancestry inference microhaplotype file generated from the step 3 (screen-MHs-by-SCGWAS) of the major pipeline.
```bash
MLE_for_MHs_genotypes -p <plink_prefix> -m <mhs_info_file> -o <output_file>
```
Intermediate files: `SNPs_in_all_candidate_MHs.snplist`, `SNPs_in_all_candidate_MHs.ped`, and `SNPs_in_all_candidate_MHs.map`.
Result file: `<output_file>`. It stores the phased MH genotype table inferred by haplo.stats.

#### `Prep_file_for_Infocalc`
Convert MH genotype file to Infocalc input format.
```bash
Prep_file_for_Infocalc -g <mhs_genotype_file> -s <sample_info_file> -o <output_file>
```
Result file: `<output_file>`. It is the Infocalc-formatted genotype input with translated MH genotype codes.

#### `Get_MHlist_with_max_In`
Keep one locus with max In for duplicated index SNP groups based on the output of Infocalc. Before this step, you should run Infocalc (https://rosenberglab.stanford.edu/infocalc.html) using the recoded MH-genotype file and get the output.
```bash
Get_MHlist_with_max_In -i <in_file> -o <output_file>
```
Result file: `<output_file>`. A MH-list contains the selected MHs with max In.

#### `Filter_MH_genotype_file_based_MHlist`
Filter an MH genotype file by a given MH list.
```bash
Filter_MH_genotype_file_based_MHlist -g <mhs_genotype_file> -m <mh_list_file> -o <output_file>
```
Result file: `<output_file>`. It is the filtered MH genotype file containing only loci from the provided MH list. As for <mhs_info_file>, the redundant MHs should be deleted artificially for subsequent steps.

#### `LOOCV_of_MHs`
Run leave-one-out cross-validation for MH-based Naive Bayes model for population classification. The input <mhs_genotype_file> is from the last step.
```bash
LOOCV_of_MHs -g <mhs_genotype_file> -s <sample_info_file> -o <output_prefix>
```
Intermediate files: `<output_prefix>_yPred.txt` and `<output_prefix>_yTrue.txt`.
Result file: `<output_prefix>_accuracy.txt`, `<output_prefix>_confusion_matrix.xlsx`, `<output_prefix>_evaluation_metrics.xlsx`, and `<output_prefix>_confusion_matrix.tiff`. They provide model accuracy, confusion matrix table, evaluation metrics, and confusion matrix plot.

#### `Naive_Bayes_random_kept_MHs`
Simulates keeping only part of the MH panel: for each possible number of retained markers, it draws that many MHs at random and runs times independent subsets (default **10000** via `-t` / `--times`). Each run uses Naive Bayes with LOOCV for ancestry prediction, so it can be observed how the panel performs as marker count changes.
```bash
Naive_Bayes_random_kept_MHs -g <mhs_genotype_file> -s <sample_info_file> -r <seed> -t <times> -d <output_dir>
```
Result file: For each `<output_dir>/kept_<n>` directory, `ACC.txt`, `AUC.txt`, `populationPredAll.txt`, and `MHs_kept_IF_all.txt`. They record simulation accuracy/AUC, per-sample predictions, and MH keep/drop indicators.

#### `Validation_of_MHs_using_test_set`
Validate MH model performance on an test set.
```bash
Validation_of_MHs_using_test_set -g <train_genotype_file> -t <test_genotype_file> -s <sample_info_file> -o <output_prefix>
```
Intermediate files: `<output_prefix>_yPred.txt` and `<output_prefix>_yTrue.txt`.
Result file: `<output_prefix>_accuracy.txt`, `<output_prefix>_confusion_matrix.xlsx`, `<output_prefix>_evaluation_metrics.xlsx`, and `<output_prefix>_confusion_matrix.tiff`. They summarize test-set classification performance.

#### `Draw_manhattan_qq_of_SCGWAS`
Draw Manhattan and QQ plots from SCGWAS result files.
```bash
Draw_manhattan_qq_of_SCGWAS -m <mhs_info_file> -p <p_value_threshold> -d <scgwas_dir> -o <output_prefix>
```
Intermediate files: `SCGWAS_<POP>_merge.assoc.logistic` (one per population).
Result file: `<output_prefix>_manhattan_QQ_plots.tiff`. This is a patchwork plot containing the Manhattan plots of the SCGWAS results for each population and the QQ plot of all SCGWAS results.

#### `Draw_PCA_of_MHs`
Draw 3D PCA plot of samples using MH genotypes.
```bash
Draw_PCA_of_MHs -g <mhs_genotype_file> -s <sample_info_file> -o <output_prefix>
```
Intermediate files: `<output_prefix>_MHs_SNPs_genotypes_for_adegenet_PCA.txt`.
Result file: `<output_prefix>_PCA_3dplot.tiff`. This is the 3D PCA plot.

#### `AUC_ACC_decline_with_n_MHs`
Plot AUC and ACC decline versus number of kept MHs from random-kept simulations and full-MH LOOCV result.
```bash
AUC_ACC_decline_with_n_MHs -d <random_kept_result_dir> -i <loocv_all_mhs_output_prefix> -o <output_prefix>
```
Result file: `<output_prefix>_AUC_ACC_decline_with_n_MHs.tiff` and `<output_prefix>_AUC_median_by_n_MHs.xlsx`. They provide the decline plot and the median/reduction summary table.

#### `Stat_MHs_indexSNPs_info`
Summarize index SNP statistics.
```bash
Stat_MHs_indexSNPs_info -m <mhs_info_file> -s <sample_info_file> -p <plink_prefix> -o <output_prefix>
```
Intermediate files: `<output_prefix>_indexSNPs.snplist` and PLINK recode file `<output_prefix>_indexSNPs.raw` (plus PLINK logs).
Result file: `<output_prefix>_indexSNPs_info.xlsx`. It reports index SNP allele, population frequencies, and multiple-regression P values.

#### `Organize_MHs_infomation_table`
Organize MH basic information into one summary table.
```bash
Organize_MHs_infomation_table -m <mhs_info_file> -p <plink_prefix> -o <output_prefix>
```
Intermediate files: `<output_prefix>_all_SNPs.snplist` and PLINK subset files such as `<output_prefix>_all_SNPs.bed/.bim/.fam` (plus PLINK logs).
Result file: `<output_prefix>_organized_MHs_info.xlsx`. It reports integrated MH information including index SNPs, Chr, MBP, constitutional unit, N SNPs, length and Ae values.

#### `Stat_allele_freq_of_MHs`
Calculate allele frequency table for each MH in each population.
```bash
Stat_allele_freq_of_MHs -g <mhs_genotype_file> -s <sample_info_file> -o <output_prefix>
```
Result file: `<output_prefix>_alleleFreqTable.xlsx`. It stores allele frequency of each MH allele by population.

#### `Calculate_Ae_of_MHs`
Calculate Ae values for MHs across populations.
```bash
Calculate_Ae_of_MHs -g <mhs_genotype_file> -s <sample_info_file> -o <output_prefix>
```
Result file: `<output_prefix>_AeTable.xlsx`. It contains Ae values of each MH across populations.

After you have run these tools in sequence as your workflow requires, you still need to look up cytoband annotations for the MHs example using [SNPnexus](https://www.snp-nexus.org/)—then format tables and assemble multi-panel figures. Together with the outputs of tools, that yields figures and tables similar to those in our article.

## Notes

1. Ensure PLINK 1.9 is properly installed and accessible from command line
2. Ensure R 4.3+ is properly installed and accessible from command line
3. Ensure R package 'haplo.stats' is installed
4. Input file paths must exist and have correct format
5. Temporary files will be generated in current directory, recommend running in dedicated working directory
6. Pay attention to memory usage when processing large files
