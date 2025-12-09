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
                    if major_version < 3 or (major_version == 3 and minor_version < 6):
                        print(f"Error: R version 3.6 or higher required, current version: {part}")
                        sys.exit(1)
                    break
            else:
                print("Warning: Could not parse R version, please ensure R 3.6+ is installed")
        else:
            print("Warning: Could not determine R version, please ensure R 3.6+ is installed")
            
        print(f"R version check passed: {version_output}")
        
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: R not installed or not accessible")
        print("Please install R 3.6 or higher and ensure it's in your PATH")
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
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.7",
    install_requires=[
        "pandas>=1.0.0",
        "openpyxl>=2.6.0",
    ],
    # Installation checks
    setup_requires=[
        "pandas>=1.0.0",
        "openpyxl>=2.6.0",
    ],
    entry_points={
        "console_scripts": [
            "snp-filter-based-pop=cli:snp_filter_main",
            "prep-MHs=cli:prep_mhs_main", 
            "screen-MHs-by-SCGWAS=cli:screen_mhs_main",
        ],
    },
    py_modules=["cli"],
    keywords="bioinformatics, genetics, microhaplotype, gwas, plink",
    project_urls={
        "Bug Reports": "https://github.com/your-repo/mh-analysis/issues",
        "Source": "https://github.com/your-repo/mh-analysis",
    },
)