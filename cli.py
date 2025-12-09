#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MHfinder
"""

import sys
import argparse
from core.snp_filter import SNPFilter
from core.prep_MHs import MicrohaplotypeGenerator
from core.scgwas_screen import SCGWASScreener

def snp_filter_main():
    """SNP filtering command line interface"""
    parser = argparse.ArgumentParser(description='Filter SNPs based on MAF in population')
    parser.add_argument('--input-prefix', '-i', type=str, required=True,
                       help='PLINK binary file prefix')
    parser.add_argument('--sample-info', '-s', type=str, required=True,
                       help='Sample information file')
    parser.add_argument('--maf-threshold', '-t', type=float, required=True,
                       help='MAF threshold')
    parser.add_argument('--filter-mode', '-m', type=str, required=True,
                       choices=['all', 'intersection', 'union'],
                       help='Filtering mode: all, intersection, union')
    parser.add_argument('--output-prefix', '-o', type=str, required=True,
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    try:
        snp_filter = SNPFilter()
        result = snp_filter.filter_snps_by_maf(
            input_prefix=args.input_prefix,
            sample_info=args.sample_info,
            maf_threshold=args.maf_threshold,
            filter_mode=args.filter_mode,
            output_prefix=args.output_prefix
        )
        print(f"SNP filtering completed: {result}")
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


def prep_mhs_main():
    """Microhaplotype preparation command line interface"""
    parser = argparse.ArgumentParser(description='Generate microhaplotypes based LD calculation result')
    parser.add_argument('--input-prefix', '-i', type=str, required=True,
                       help='PLINK binary file prefix')
    parser.add_argument('--snplist', '-s', type=str, required=True,
                       help='Snplist file from previous step')
    parser.add_argument('--length-threshold', '-l', type=int, required=True,
                       help='Length threshold (distance threshold between two endpoint SNPs of microhaplotype, unit bp)')
    parser.add_argument('--min-site-number', '-n', type=int, required=True,
                       help='Minimum number of SNPs in microhaplotype')
    parser.add_argument('--r2-threshold', '-r', type=float, required=True,
                       help='Strong linkage threshold (r2 between any two SNPs in microhaplotype must be below this value)')
    parser.add_argument('--output', '-o', type=str, required=True,
                       help='Output file name')
    
    args = parser.parse_args()
    
    try:
        mh_generator = MicrohaplotypeGenerator()
        result = mh_generator.prepare_microhaplotypes(
            input_prefix=args.input_prefix,
            snplist=args.snplist,
            length_threshold=args.length_threshold,
            min_site_number=args.min_site_number,
            r2_threshold=args.r2_threshold,
            output=args.output
        )
        print(f"Microhaplotype preparation completed: {result}")
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


def screen_mhs_main():
    """Stepwise conditional GWAS screening command line interface"""
    parser = argparse.ArgumentParser(description='Screen microhaplotypes based on conditional GWAS')
    parser.add_argument('--input-prefix', '-i', type=str, required=True,
                       help='PLINK binary file prefix')
    parser.add_argument('--sample-info', '-s', type=str, required=True,
                       help='Sample information file')
    parser.add_argument('--mhs-file', '-m', type=str, required=True,
                       help='Microhaplotype file from previous step')
    parser.add_argument('--p-threshold', '-p', type=float, required=True,
                       help='GWAS p-value threshold')
    parser.add_argument('--mhs-diversity-threshold', '-d', type=float, required=False,
                       help='Diversity threshold of microhaplotype (optional)')
    parser.add_argument('--Ae-threshold', '-a', type=float, required=False,
                       help='Effective number of alleles threshold (optional)')
    parser.add_argument('--output', '-o', type=str, required=True,
                       help='Output file name')
    parser.add_argument('--loop-times', '-n', type=int,
                       help='Maximum GWAS loop times (optional, default unlimited)')
    
    args = parser.parse_args()
    
    try:
        scgwas_screener = SCGWASScreener()
        result = scgwas_screener.screen_mhs_by_scgwas(
            input_prefix=args.input_prefix,
            sample_info=args.sample_info,
            mhs_file=args.mhs_file,
            p_threshold=args.p_threshold,
            mhs_diversity_threshold=args.mhs_diversity_threshold,
            Ae_threshold=args.Ae_threshold,
            output=args.output,
            loop_times=args.loop_times
        )
        print(f"Stepwise conditional GWAS screening completed: {result}")
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)
