#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MHfinder
"""

import subprocess
import sys
import platform
from setuptools import setup, find_packages

def check_plink_version():
    """Check PLINK version during installation"""
    try:
        result = subprocess.run(['plink', '--version'], 
                              capture_output=True, text=True, check=True)
        version_output = result.stdout.strip()
        
        if not version_output.startswith('PLINK v1.9'):
            print(f"Error: PLINK v1.9 required, current version: {version_output}")
            sys.exit(1)
            
        print(f"PLINK version check passed: {version_output}")
        
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: PLINK not installed or not accessible")
        sys.exit(1)

def check_r_version():
    """Check R version during installation"""
    try:
        result = subprocess.run(['Rscript', '--version'], 
                              capture_output=True, text=True, check=True)
        version_output = result.stdout.strip()
        
        # Extract R version number
        version_line = [line for line in version_output.split('\n') if 'R version' in line]
        if version_line:
            version_str = version_line[0]
            # Extract version number (e.g., "3.6.0" from "R version 3.6.0")
            version_parts = version_str.split()
            for part in version_parts:
                if part.count('.') >= 2 and part.replace('.', '').isdigit():
                    major_version = int(part.split('.')[0])
                    minor_version = int(part.split('.')[1])
                    if major_version < 4 or (major_version == 4 and minor_version < 3):
                        print(f"Error: R version 4.3 or higher required, current version: {part}")
                        sys.exit(1)
                    break
            else:
                print("Warning: Could not parse R version, please ensure R 4.3+ is installed")
        else:
            print("Warning: Could not determine R version, please ensure R 4.3+ is installed")
            
        print(f"R version check passed: {version_output}")
        
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: R not installed or not accessible")
        print("Please install R 4.3 or higher and ensure it's in your PATH")
        sys.exit(1)

def check_r_package():
    """Check if haplo.stats package is installed in R"""
    try:
        result = subprocess.run(['Rscript', '-e', '"library(haplo.stats)"'], 
                              capture_output=True, text=True, check=True)
        print("R package haplo.stats check passed")
    except subprocess.CalledProcessError:
        print("Error: R package 'haplo.stats' not found")
        print("Please install it in R using: install.packages('haplo.stats')")
        sys.exit(1)

# Check dependencies during installation
if 'install' in sys.argv or 'develop' in sys.argv:
    check_plink_version()
    check_r_version()
    check_r_package()

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="MHfinder",
    version="1.0.0",
    author="Microhaplotype Analysis Team",
    description="Microhaplotype analysis module integrating SNP filtering based on MAF in population, microhaplotype preparation and conditional GWAS screening",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    install_requires=[
        "pandas>=2.0.0",
        "openpyxl>=3.1.0",
    ],
    # Installation checks
    setup_requires=[
        "pandas>=2.0.0",
        "openpyxl>=3.1.0",
    ],
    entry_points={
        "console_scripts": [
            "snp-filter-based-pop=cli:snp_filter_main",
            "prep-MHs=cli:prep_mhs_main", 
            "screen-MHs-by-SCGWAS=cli:screen_mhs_main",
            "MLE_for_MHs_genotypes=extra.cli_entrypoints:mle_mh_genotypes_cli",
            "Draw_manhattan_qq_of_SCGWAS=extra.cli_entrypoints:plot_scgwas_manhattan_qq_cli",
            "Draw_PCA_of_MHs=extra.cli_entrypoints:plot_mh_pca_cli",
            "LOOCV_of_MHs=extra.cli_entrypoints:loocv_mhs_cli",
            "Validation_of_MHs_using_test_set=extra.cli_entrypoints:validate_mhs_testset_cli",
            "Organize_MHs_infomation_table=extra.cli_entrypoints:organize_mh_info_table_cli",
            "Stat_allele_freq_of_MHs=extra.cli_entrypoints:stat_mh_allele_freq_cli",
            "Stat_MHs_indexSNPs_info=extra.cli_entrypoints:stat_mh_indexsnps_cli",
            "Calculate_Ae_of_MHs=extra.cli_entrypoints:calc_mh_ae_cli",
            "Naive_Bayes_random_kept_MHs=extra.cli_entrypoints:naivebayes_random_kept_mhs_cli",
            "Filter_MH_genotype_file_based_MHlist=extra.cli_entrypoints:filter_mh_genotype_by_list_cli",
            "Get_MHlist_with_max_In=extra.cli_entrypoints:get_mhlist_max_in_cli",
            "Prep_file_for_Infocalc=extra.cli_entrypoints:prep_infocalc_input_cli",
            "AUC_ACC_decline_with_n_MHs=extra.cli_entrypoints:auc_acc_decline_with_n_mhs_cli",
        ],
    },
    py_modules=["cli"],
    keywords="bioinformatics, genetics, microhaplotype, gwas, plink",
    project_urls={
        "Bug Reports": "https://github.com/your-repo/mh-analysis/issues",
        "Source": "https://github.com/your-repo/mh-analysis",
    },
)
