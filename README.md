# MHfinder - Microhaplotype Analysis Module

## Project Overview

This is a Python module for microhaplotype analysis, integrating three core functions:
1. **Population-based MAF SNP filtering** - Filter SNPs based on Minor Allele Frequency (MAF)
2. **Microhaplotype generation** - Generate microhaplotypes based on filtered SNPs
3. **Stepwise conditional GWAS-based microhaplotype screening** - Screen candidate microhaplotypes using stepwise conditional GWAS method

## Project Structure

```
buildModule/
├── core/                    # Core functionality modules
│   ├── __init__.py         # Module initialization
│   ├── utils.py            # General utility functions
│   ├── snp_filter.py       # SNP filtering functionality
│   ├── prep_MHs.py         # Microhaplotype generation functionality
│   └── scgwas_screen.py     # Stepwise conditional GWAS screening functionality
├── cli.py                  # Command line interface
├── README.md               # Project documentation
├── requirements.txt        # Dependency package list
├── setup.py               # Installation configuration
└── test_module.py         # Test module
```

## Requirements

- **PLINK 1.9**: PLINK 1.9 version must be installed, check version by running `plink --version`
- **R >= 3.6**: R version 3.6 or higher must be installed
- **Python >= 3.7**: Python version 3.7 or higher required

### R Package Requirements
- **haplo.stats**: Required R package for haplotype analysis
  ```r
  install.packages("haplo.stats")
  ```

### Python Package Requirements
- **pandas >= 1.0.0**: Data manipulation and analysis
- **openpyxl >= 2.6.0**: Excel file support

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

### 2. Install R (>= 3.6)
```bash
# Windows: Download from https://cran.r-project.org/bin/windows/base/
# macOS: Download from https://cran.r-project.org/bin/macosx/
# Linux: Use package manager (e.g., sudo apt-get install r-base)

# Verify installation
Rscript --version
# Should display: R version 3.6.x or higher
```

### 3. Install Required R Package
```r
# Start R and install the required package
R
install.packages("haplo.stats")
```

### 4. Install Python Module
```bash
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
- `mhs_diversity_threshold`: minimum Diversity value of microhaplotype
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
    --input-prefix data/genotype \
    --sample-info data/sample_info.csv \
    --maf-threshold 0.1 \
    --filter-mode union \
    --output-prefix filtered_snps

# 2. Generate microhaplotypes
prep_MHs \
    --input-prefix data/genotype \
    --snplist filtered_snps.snplist \
    --length-threshold 150 \
    --min-site-number 3 \
    --r2-threshold 0.8 \
    --output microhaplotypes.txt

# 3. Screen microhaplotypes based on stepwise conditional GWAS
screen-MHs-by-SCGWAS \
    --input-prefix data/genotype \
    --sample-info data/sample_info.csv \
    --mhs-file microhaplotypes.txt \
    --p-threshold 5e-8 \
    --mhs-diversity-threshold 0.6 \
    --Ae-threshold 3.0 \
    --loop-times 5 \
    --output final_microhaplotypes.txt
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

## Notes

1. Ensure PLINK 1.9 is properly installed and accessible from command line
2. Ensure R 3.6+ is properly installed and accessible from command line
3. Ensure R package 'haplo.stats' is installed
4. Input file paths must exist and have correct format
5. Temporary files will be generated in current directory, recommend running in dedicated working directory
6. Pay attention to memory usage when processing large files
