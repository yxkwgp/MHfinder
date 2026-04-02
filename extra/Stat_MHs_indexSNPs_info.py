import argparse
import os
import pandas as pd
import numpy as np
import statsmodels.api as sm

parser = argparse.ArgumentParser(description = 'Stat basic information of index SNPs in MHs')
parser.add_argument('--mhs-info-file', '-m', type=str, required=True,
                       help='MHs information file')
parser.add_argument('--sample-info-file', '-s', type=str, required=True,
                       help='Sample information file')
parser.add_argument('--plink-prefix', '-p', type=str, required=True,
                       help='Plink binary file prefix')
parser.add_argument('--output-prefix', '-o', type=str, required=True,
                       help='Output filename prefix')
args = parser.parse_args()
mhs_info_file = args.mhs_info_file
sample_info_file = args.sample_info_file
plink_prefix = args.plink_prefix
output_prefix = args.output_prefix

# Get index SNPs of each population
i = 1
pop_indexSNPs = dict()
infile1 = open(mhs_info_file, 'r')
for line in infile1:
    if i % 3 == 1:
        line = line.strip().lstrip('>').split('_')
        pop = line[0]
        if pop not in pop_indexSNPs.keys():
            pop_indexSNPs[pop] = []
    elif i % 3 == 2:
        line = line.strip().split('\t')
        indexSNP = line[0].lstrip('*')
        pop_indexSNPs[pop].append(indexSNP)
    i += 1
infile1.close()

# Extract index SNPs from PLINK file
outfile1 = open(f'{output_prefix}_indexSNPs.snplist', 'w')
for pop in pop_indexSNPs.keys():
    outfile1.write('\n'.join(pop_indexSNPs[pop]) + '\n')
outfile1.close()
cmd = f'plink --bfile {plink_prefix} --extract {output_prefix}_indexSNPs.snplist --recodeA --out {output_prefix}_indexSNPs'
os.system(cmd)

# Read in .raw, sample info and get information
raw_file_name = f'{output_prefix}_indexSNPs.raw'
raw_df = pd.read_csv(raw_file_name, sep = ' ')
sample_info_df = pd.read_csv(sample_info_file, sep = '\t')
sample_info_df.columns = ['FID', 'IID', 'population']

## Allele
SNPs_colnames = raw_df.columns[6:]
SNPs_alleles_dict = dict()
for SNP_colname in SNPs_colnames:
    SNP_colname = SNP_colname.split('_')
    SNPs_alleles_dict[SNP_colname[0]] = SNP_colname[1]
SNPs_alleles_df = pd.DataFrame(data = SNPs_alleles_dict.items(), columns = ['SNP', 'Allele'])

## Allele frequency
raw_df_withInfo = pd.merge(raw_df, sample_info_df[['IID', 'population']], on = 'IID', how = 'left').drop(columns = ['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])
raw_df_withInfo = raw_df_withInfo[['IID', 'population'] + list(SNPs_colnames)]
SNPs_colnames_new = [new_colname.split('_')[0] for new_colname in raw_df_withInfo.columns[2:]]
raw_df_withInfo.columns = ['IID', 'population'] + SNPs_colnames_new
raw_df_withInfo_group = raw_df_withInfo.groupby('population')
popOrder = raw_df_withInfo['population'].drop_duplicates().sort_values().tolist()
allele_freq_table = np.zeros(shape = (len(raw_df_withInfo.columns) - 2, len(popOrder)))
allele_freq_table = pd.DataFrame(data = allele_freq_table, index = SNPs_colnames_new, columns = popOrder)
for pop in popOrder:
    raw_df_part = raw_df_withInfo_group.get_group(pop).drop(columns = ['IID', 'population'])
    for i in range(raw_df_part.shape[1]):
        allele_count = raw_df_part.iloc[:, i].sum()
        allele_freq = allele_count / (raw_df_part.shape[0] * 2)
        allele_freq_table.loc[SNPs_colnames_new[i], pop] = allele_freq
allele_freq_table = allele_freq_table.reset_index()
allele_freq_table = allele_freq_table.rename(columns={allele_freq_table.columns[0]: 'SNP'})

## Multiple regression P
all_multiple_P_list = []
for pop in popOrder:
    indexSNPs = pop_indexSNPs[pop]
    raw_df_vertical_part = raw_df_withInfo[['IID', 'population'] + indexSNPs].copy()
    raw_df_vertical_part['IF'] = (raw_df_vertical_part['population'] == pop).astype(int)
    y = raw_df_vertical_part['IF']
    X = raw_df_vertical_part[indexSNPs]
    X = sm.add_constant(X)
    model = sm.GLM(y, X, family=sm.families.Binomial())
    fit_result = model.fit()
    indexSNPs_multiP = fit_result.summary2().tables[1].loc[indexSNPs, 'P>|z|'].reset_index()
    indexSNPs_multiP.columns = ['SNP', 'Multiple regression P']
    all_multiple_P_list.append(indexSNPs_multiP)
all_multiple_P_df = pd.concat(all_multiple_P_list, ignore_index = True)

## Merge all information
result_df = pd.merge(all_multiple_P_df, SNPs_alleles_df, on = 'SNP', how = 'left')
result_df = pd.merge(result_df, allele_freq_table, on = 'SNP', how = 'left')
result_df.to_excel(f'{output_prefix}_indexSNPs_info.xlsx', index = False)
