#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
General utility functions module
"""

import os
import subprocess
import pandas as pd

def read_table(filename):
    """Read xlsx, csv, tsv, or newline-separated txt files into pandas DataFrame"""
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    
    try:
        if filename.endswith('.xlsx'):
            df = pd.read_excel(filename)
        elif filename.endswith('.csv'):
            df = pd.read_csv(filename, header=0)
        elif filename.endswith('.tsv') or filename.endswith('.txt'):
            df = pd.read_table(filename, header=0)
        else:
            raise ValueError(f'Unsupported file format: {filename}, file extension must be xlsx, csv, tsv or txt')
        return df
    except Exception as e:
        raise RuntimeError(f"Failed to read file {filename}: {str(e)}")


def run_plink_command(cmd):
    """Run PLINK command and check result"""
    print(f"Executing command: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"PLINK command execution failed: {result.stderr}")
    return result.stdout


def run_rscript_command(r_script, *args):
    """Run Rscript command and check result"""
    cmd = ['Rscript', r_script] + list(args)
    print(f"Executing R command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Rscript command execution failed: {result.stderr}")
    return result.stdout 