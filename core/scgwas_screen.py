#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stepwise conditional GWAS screening functionality module
"""

import os
import pandas as pd
import numpy as np
from gc import collect
from .utils import read_table, run_plink_command, run_rscript_command

class SCGWASScreener:
    """Stepwise conditional GWAS screening class"""
    
    def _find_mhs(self, index_snp, mhs_list):
        """Find all MHs containing indexSNP"""
        selected_mhs = []
        for element in mhs_list:
            if index_snp in element:
                selected_mhs.append(list(element))
        return selected_mhs

    def _GTbyHS(self, geno_df):
        """Convert genotype to haplotype using Rscript"""
        input_file = 'temp_haplotype_input.tsv'
        output_file = 'temp_haplotype_output.tsv'
        geno_df.to_csv(input_file, sep='\t', index=False)
        r_script_path = os.path.join(os.path.dirname(__file__), '..', 'MH_phasing.R')
        run_rscript_command(r_script_path, input_file, output_file)
        result_df = pd.read_csv(output_file, sep='\t')
        return result_df
  
    def _check_diversity(self, mhs_list, input_prefix, mhs_diversity_threshold):
        """Check if MHs Diversity value is greater than or equal to threshold"""
        checked_d_value = []
        diversity_checked = []
        snps = set()
        for mh in mhs_list:
            snps.update(mh)
        
        cmd = f'plink --bfile {input_prefix} --snps {", ".join(snps)} --recodeA --out tem_file_for_calculate_Diversity'
        run_plink_command(cmd)
        
        infile_raw = open('tem_file_for_calculate_Diversity.raw', 'r')
        first_line = infile_raw.readline().strip().split()
        infile_raw.close()
        
        raw_colnames = first_line[0:6]
        for element in first_line[6:]:
            element = element.split('_')
            raw_colnames.append(element[0])
        
        geno_number = pd.read_csv('tem_file_for_calculate_Diversity.raw', sep=' ')
        geno_number.columns = raw_colnames
        
        for mh in mhs_list:
            upstairs = 0
            downstairs = geno_number.shape[0] * (geno_number.shape[0] - 1) / 2
            
            for m in range(geno_number.shape[0] - 1):
                for n in range(m + 1, geno_number.shape[0]):
                    count_snp = 0
                    for snp in mh:
                        if geno_number[snp].iloc[m] != geno_number[snp].iloc[n]:
                            count_snp += 1
                    if count_snp >= 2:
                        upstairs += 1
            
            diversity = upstairs / downstairs
            if diversity >= mhs_diversity_threshold:
                checked_d_value.append(mh)
                diversity_checked.append(diversity)
        
        del geno_number
        if os.path.exists('tem_file_for_calculate_Diversity.raw'):
            os.remove('tem_file_for_calculate_Diversity.raw')
        collect()
        
        return checked_d_value, diversity_checked

    def _check_Ae(self, mhs_list, input_prefix, sample_info, Ae_threshold):
        """Check if MHs Effective number of alleles is greater than or equal to threshold"""
        checked_Ae_value = []
        Ae_checked = []
        snps = set()
        for mh in mhs_list:
            snps.update(mh)
        cmd = f'plink --bfile {input_prefix} --snps {", ".join(snps)} --recode --out tem_file_for_calculate_Ae'
        run_plink_command(cmd)
        dfPed = pd.read_csv('tem_file_for_calculate_Ae.ped', sep=' ', header=None)
        dfMap = pd.read_csv('tem_file_for_calculate_Ae.map', header=None, sep = '\t')
        dfMap['a1'] = dfMap[1] + '.a1'
        dfMap['a2'] = dfMap[1] + '.a2'
        dfPed.columns = ['FID', 'IID', 'mother', 'father', 'sex', 'pheno'] + list(dfMap[['a1', 'a2']].to_numpy().flatten(order='C'))
        for mh in mhs_list:
            label = [f'{snp}.a{allele}' for snp in mh for allele in [1, 2]]
            dfGeno = dfPed[label]
            MH_colname = ['tmp.A1', 'tmp.A2']
            dfPhased = self._GTbyHS(dfGeno)
            dfPhased = dfPhased.reset_index(drop = True)
            sample_info.columns = ['IID', 'population']
            dfPhased = pd.concat([dfPhased, dfPed[['IID']]], axis=1)
            dfPhased = dfPhased.merge(sample_info, on = 'IID', how = 'inner')
            ae_values = []
            ae_values_rem = []
            for population, dfPhased_sub in dfPhased.groupby('population'):
                genotypes = dfPhased_sub.iloc[:,:2].values.flatten()
                genotypes_counts = pd.Series(genotypes).value_counts()
                total_counts = len(genotypes)
                genotype_frequencies = genotypes_counts / total_counts
                sum_squared_frequencies = (genotype_frequencies ** 2).sum()
                ae_value = 1 / sum_squared_frequencies
                ae_value_rem = f'Ae({population})={ae_value}'
                ae_values.append(ae_value)
                ae_values_rem.append(ae_value_rem)
            if np.mean(ae_values) >= Ae_threshold:
                checked_Ae_value.append(mh)
                Ae_checked.append(', '.join(ae_values_rem))
        return checked_Ae_value, Ae_checked

    def screen_mhs_by_scgwas(self, input_prefix, sample_info, mhs_file, p_threshold, 
                           output, mhs_diversity_threshold=None, Ae_threshold=None, loop_times=None):
        """
        Screen microhaplotypes based on stepwise conditional GWAS
        
        Parameters:
            input_prefix (str): PLINK binary file prefix
            sample_info (str): Sample information file
            mhs_file (str): Microhaplotype file from previous step
            p_threshold (float): GWAS p-value threshold
            mhs_diversity_threshold (float, optional): Diversity threshold of microhaplotype
            Ae_threshold (float, optional): Effective number of alleles threshold
            output (str): Output file name
            loop_times (int, optional): Maximum GWAS loop times, default unlimited
            
        Returns:
            str: Screened microhaplotype file path
        """
        # Parameter validation
        if not os.path.exists(f"{input_prefix}.bed"):
            raise FileNotFoundError(f"PLINK file not found: {input_prefix}.bed")
        if not os.path.exists(mhs_file):
            raise FileNotFoundError(f"Microhaplotype file not found: {mhs_file}")
        if p_threshold <= 0 or p_threshold > 1:
            raise ValueError("p_threshold must be between 0 and 1")
        if mhs_diversity_threshold is not None and (mhs_diversity_threshold < 0 or mhs_diversity_threshold > 1):
            raise ValueError("mhs_diversity_threshold must be between 0 and 1")
        if Ae_threshold is not None and Ae_threshold < 1:
            raise ValueError("Ae_threshold must be greater than or equal to 1")
        if loop_times is not None and loop_times <= 0:
            raise ValueError("loop_times must be greater than 0")
        if loop_times is None:
            loop_times = float('inf')
        
        print(f"Starting stepwise conditional GWAS-based microhaplotype screening, p threshold: {p_threshold}, diversity threshold: {mhs_diversity_threshold}, Ae threshold: {Ae_threshold}")
        
        # Read all MHs and generate SNP list
        mhs_snp_set = set()
        mhs_list = []
        mhs_infile = open(mhs_file, 'r')
        mhs_infile.readline()
        for line in mhs_infile:
            line = line.strip().split('|')
            snp_list = line[1].split('-')
            mhs_snp_set.update(snp_list)
            mhs_list.append(tuple(snp_list))
        mhs_infile.close()
        mhs_list = tuple(mhs_list)
        
        mhs_snp_outfile = open('tmp.snplist', 'w')
        for snp in mhs_snp_set:
            mhs_snp_outfile.write(snp + '\n')
        mhs_snp_outfile.close()
        
        # Execute recoding
        sample_info_df = read_table(sample_info)
        sample_genotype_col = sample_info_df.columns[2]
        sample_genotype_all = sorted(sample_info_df[sample_genotype_col].unique().tolist())
        
        sample_recoded_df = pd.DataFrame()
        sample_recoded_df['FID'] = sample_info_df[sample_info_df.columns[0]]
        sample_recoded_df['IID'] = sample_info_df[sample_info_df.columns[1]]
        
        for genotype in sample_genotype_all:
            sample_recoded_df[genotype] = (sample_info_df[sample_genotype_col] == genotype).astype(int) + 1
        
        sample_recoded_df.to_csv('tmp.sample', sep='\t', index=False)
        
        # Multiple SCGWAS loops
        output_file = open(output, 'w')
        
        for genotype in sample_genotype_all:
            print(f"Processing genotype: {genotype}")
            
            # First GWAS
            cmd = f'plink --bfile {input_prefix} --extract tmp.snplist --pheno tmp.sample --pheno-name {genotype} --logistic hide-covar --ci 0.95 --allow-no-sex --out SCGWAS_{genotype}_1'
            run_plink_command(cmd)
            
            index_snps = []
            exclud_snps = set()
            i = 1
            loop_clock = loop_times
            
            while loop_clock > 0:
                # Process PLINK result file line by line with Python, filter lines with p-value less than threshold
                infile = open(f'SCGWAS_{genotype}_{i}.assoc.logistic', 'r')
                infile.readline()
                outfile = open(f'SCGWAS_{genotype}_{i}.filter.assoc.logistic', 'w')
                for line in infile:
                    line = line.strip().split()
                    if line[11] != 'NA' and float(line[11]) < p_threshold:
                        outfile.write(f'{line[1]} {line[0]} {line[2]} {line[11]}\n')
                infile.close()
                outfile.close()
                data_gwas = pd.read_csv(f'SCGWAS_{genotype}_{i}.filter.assoc.logistic', sep=' ', header=None, names=['SNP', 'CHR', 'BP', 'P'], dtype={'SNP': str, 'CHR': str, 'BP': str, 'P': float})
                if data_gwas.shape[0] == 0:
                    print(f'{genotype} Loop {i}: SCGWAS_{genotype}_{i}.filter.assoc.logistic is empty, the process for {genotype} is finished')
                    break
                
                data_gwas = data_gwas.sort_values(by='P')
                
                for j in range(data_gwas.shape[0]):
                    index_snp = data_gwas.iloc[j, 0]
                    index_snp_chr = data_gwas.iloc[j, 1]
                    print(f'{genotype} Loop {i}: {index_snp} is selected')
                    
                    selected_mhs = self._find_mhs(index_snp, mhs_list)
                    print(f'{genotype} Loop {i}: {index_snp} in MHs: {selected_mhs}')
                    
                    if mhs_diversity_threshold is not None and Ae_threshold is not None:
                        mhs_diversity_checked, diversity_checked = self._check_diversity(selected_mhs, input_prefix, mhs_diversity_threshold)
                        if len(mhs_diversity_checked) == 0:
                            mhs_checked = mhs_diversity_checked
                        else:
                            mhs_checked, Ae_value_checked = self._check_Ae(mhs_diversity_checked, input_prefix, sample_info_df.iloc[:, 1:3], Ae_threshold)
                            diversity_checked = [diversity_checked[mhs_diversity_checked.index(MH)] for MH in mhs_checked]
                    elif mhs_diversity_threshold is not None:
                        mhs_checked, diversity_checked = self._check_diversity(selected_mhs, input_prefix, mhs_diversity_threshold)
                    elif Ae_threshold is not None:
                        mhs_checked, Ae_value_checked = self._check_Ae(selected_mhs, input_prefix, sample_info_df.iloc[:, 1:3], Ae_threshold)
                    else:
                        mhs_checked = selected_mhs
                    
                    if len(mhs_checked) == 0:
                        print(f'{genotype} Loop {i}: no MHs of {index_snp} passing Diversity or Aecheck, {index_snp} is excluded')
                        exclud_snps.add(index_snp)
                        continue
                    else:
                        print(f'{genotype} Loop {i}: {mhs_checked} of {index_snp} passing Diversity or Ae check, add {index_snp} to corvariant and process next loop')
                        index_snps.append(index_snp)
                        
                        for k in range(len(mhs_checked)):
                            dict_key = f'{genotype}_{i}_{k + 1}'
                            if mhs_diversity_threshold is not None and Ae_threshold is not None:
                                output_file.write('>' + dict_key + '\n*' + index_snp + '\t' + index_snp_chr + '\tDiversity=' + str(round(diversity_checked[k], 3)) + '\t' + Ae_value_checked[k] + '\n' + '\t'.join(mhs_checked[k]) + '\n')
                            elif mhs_diversity_threshold is not None:
                                output_file.write('>' + dict_key + '\n*' + index_snp + '\t' + index_snp_chr + '\tDiversity=' + str(round(diversity_checked[k], 3)) + '\n' + '\t'.join(mhs_checked[k]) + '\n')
                            elif Ae_threshold is not None:
                                output_file.write('>' + dict_key + '\n*' + index_snp + '\t' + index_snp_chr + '\t' + Ae_value_checked[k] + '\n' + '\t'.join(mhs_checked[k]) + '\n')
                            else:
                                output_file.write('>' + dict_key + '\n*' + index_snp + '\t' + index_snp_chr + '\n' + '\t'.join(mhs_checked[k]) + '\n')
                            exclud_snps.update(mhs_checked[k])
                        
                        del data_gwas
                        collect()
                        break
                else:
                    print(f'{genotype} Loop {i}: No SNPs having MHs passing Diversity or Ae check, the process for {genotype} is finished')
                    break
                
                # Write indexSNPs to covariate file, write SNPs to be excluded to snplist
                snp_file = open(f'SCGWAS_{genotype}_{i}_indexSNPs.snplist', 'w')
                for snp in index_snps:
                    snp_file.write(snp + '\n')
                snp_file.close()
                
                exclud_snp_file = open(f'SCGWAS_{genotype}_{i}_excludSNPs.snplist', 'w')
                for snp in exclud_snps:
                    exclud_snp_file.write(snp + '\n')
                exclud_snp_file.close()
                
                cmd = f'plink --bfile {input_prefix} --extract SCGWAS_{genotype}_{i}_indexSNPs.snplist --recodeA --out SCGWAS_{genotype}_{i}_cov'
                run_plink_command(cmd)
                
                covar_file = open(f'SCGWAS_{genotype}_{i}_cov.raw', 'r')
                first_line = covar_file.readline().strip().split()
                covar_name = ' '.join(first_line[6:])
                covar_file.close()
                
                cmd = f'plink --bfile {input_prefix} --extract tmp.snplist --pheno tmp.sample --pheno-name {genotype} --exclude SCGWAS_{genotype}_{i}_excludSNPs.snplist --logistic hide-covar --ci 0.95 --allow-no-sex --covar SCGWAS_{genotype}_{i}_cov.raw --covar-name {covar_name} --out SCGWAS_{genotype}_{i+1}'
                run_plink_command(cmd)
                
                loop_clock -= 1
                i += 1
        
        output_file.close()
        
        print(f"Stepwise conditional GWAS screening completed, output file: {output}")
        return output 