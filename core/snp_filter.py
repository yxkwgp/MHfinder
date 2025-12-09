#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SNP filtering functionality module
"""

import os
from .utils import read_table, run_plink_command


class SNPFilter:
    """SNP filtering class"""
    
    def filter_snps_by_maf(self, input_prefix, sample_info, maf_threshold, 
                          filter_mode, output_prefix):
        """
        Filter SNPs based on MAF in population
        
        Parameters:
            input_prefix (str): PLINK binary file prefix
            sample_info (str): Sample information file
            maf_threshold (float): MAF threshold
            filter_mode (str): Filtering mode ('all', 'intersection', 'union')
            output_prefix (str): Output file prefix
            
        Returns:
            str: Generated snplist file
        """
        # Parameter validation
        if not os.path.exists(f"{input_prefix}.bed"):
            raise FileNotFoundError(f"PLINK file not found: {input_prefix}.bed")
        if filter_mode not in ['all', 'intersection', 'union']:
            raise ValueError("filter_mode must be 'all', 'intersection' or 'union'")
        if maf_threshold <= 0 or maf_threshold >= 0.5:
            raise ValueError("maf_threshold must be between 0 and 0.5")
        
        print(f"Starting MAF-based SNP filtering, mode: {filter_mode}")
        
        # If filterMode is 'all', filter SNPs by MAF across all samples
        if filter_mode == 'all':
            cmd = f'plink --bfile {input_prefix} --maf {maf_threshold} --write-snplist --out {output_prefix}'
            run_plink_command(cmd)
            output_file = f'{output_prefix}.snplist'
        else:
            # Create files for each population
            group_file = read_table(sample_info)
            # Use the third column as population column (index 2)
            population_col = group_file.columns[2]
            populations = group_file[population_col].unique()
            print(f"Found populations: {populations}")
            
            for pop in populations:
                pop_data = group_file[group_file[population_col] == pop]
                output_data = pop_data[[group_file.columns[0], group_file.columns[1]]].copy()  # FID and IID columns
                output_data.columns = ['FID', 'IID']
                filename = f'kept_fam_{pop}.txt'
                output_data.to_csv(filename, sep=' ', header=False, index=False)
            
            # Generate and run PLINK commands for each population to get snplist based on MAF
            for pop in populations:
                cmd = f'plink --bfile {input_prefix} --keep-fam kept_fam_{pop}.txt --maf {maf_threshold} --write-snplist --out {output_prefix}_{pop}'
                run_plink_command(cmd)
            
            # Read SNP lists from all populations
            snp_sets = []
            for pop in populations:
                snplist_file = f'{output_prefix}_{pop}.snplist'
                infile = open(snplist_file, 'r')
                snps = set(line.strip() for line in infile if line.strip())
                infile.close()
                snp_sets.append(snps)
            
            # Merge SNPs list based on filterMode
            if filter_mode == 'intersection':
                # Take intersection, SNPs that satisfy MAF threshold in all populations
                if snp_sets:
                    merged_snps = snp_sets[0]
                    for snp_set in snp_sets[1:]:
                        merged_snps = merged_snps.intersection(snp_set)
                else:
                    merged_snps = set()
            elif filter_mode == 'union':
                # Take union, SNPs that satisfy MAF threshold in at least one population
                merged_snps = set()
                for snp_set in snp_sets:
                    merged_snps = merged_snps.union(snp_set)
            
            # Write merged SNP list to snplist file
            output_file = f'{output_prefix}.snplist'
            outfile = open(output_file, 'w')
            for snp in sorted(merged_snps):
                outfile.write(f'{snp}\n')
            outfile.close()
        
        print(f"SNP filtering completed, output file: {output_file}")
        return output_file 