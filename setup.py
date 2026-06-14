#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Setuptools configuration for MHfinder.

External executables such as PLINK and R are checked at run time rather than at
installation time. This keeps installation lightweight for package managers,
continuous integration, and reviewers who only want to inspect the Python API.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="MHfinder",
    version="1.0.0",
    author="Di Yu, Qian Wang, Fan Liu and colleagues",
    description=(
        "Population-aware microhaplotype discovery and stepwise conditional "
        "GWAS screening for ancestry-informative marker panel design"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    py_modules=["cli"],
    python_requires=">=3.8",
    install_requires=[
        "pandas>=2.0.0",
        "openpyxl>=3.1.0",
    ],
    extras_require={
        "analysis": [
            "numpy>=1.23",
            "scipy>=1.9",
            "scikit-learn>=1.2",
            "statsmodels>=0.13",
        ],
        "dev": [
            "pytest>=7",
        ],
    },
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
            "Organize_MHs_information_table=extra.cli_entrypoints:organize_mh_info_table_cli",
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
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    keywords=(
        "bioinformatics population-genetics microhaplotype ancestry-inference "
        "GWAS PLINK forensic-genomics marker-selection"
    ),
    project_urls={
        "Source": "https://github.com/Fun-Gene/MHfinder",
        "Bug Reports": "https://github.com/Fun-Gene/MHfinder/issues",
        "Documentation": "https://github.com/Fun-Gene/MHfinder#readme",
    },
)
