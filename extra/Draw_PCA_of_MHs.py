import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description = 'Draw PCA plot of MHs')
parser.add_argument('--mhs-genotype', '-g', type=str, required=True,
                       help='MHs genotype file')
parser.add_argument('--sample-info', '-s', type=str, required=True,
                       help='Sample information file')
parser.add_argument('--output-prefix', '-o', type=str, required=True,
                       help='Output filename prefix')
args = parser.parse_args()
mhs_genotype = args.mhs_genotype
sample_info = args.sample_info
output_prefix = args.output_prefix

def IF_split_count_consistent(pd_series):
    split_counts = pd_series.str.split('-').str.len()
    unique_counts = split_counts.unique()
    if len(unique_counts) == 1:
        return True
    else:
        return False

def Split_MHs_columns(pd_series):
    base_name = pd_series.name
    pd_series = pd_series.str.split('-', expand = True)
    n_SNPs = pd_series.shape[1]
    pd_series.columns = [f'{base_name}.{base_name}_SNP{i + 1}' for i in range(n_SNPs)]
    return pd_series

def Get_major_allele(pd_series):
    SNP_name = pd_series.name
    allele_counts = pd_series.value_counts()
    if allele_counts.size > 2:
        raise ValueError(f'{SNP_name} has more than 2 alleles')
    return allele_counts.index[0]

def Translate_genotype_to_number(pd_series, major_allele):
    new_series = pd_series.groupby(pd_series.iloc[:, 0]).apply(lambda g: (g.iloc[:, 1] != major_allele).sum())
    new_series.reset_index(drop = True, inplace = True)
    return new_series

# Reshape MH-genotype dataframe        
MHs_genotype_df = pd.read_csv(mhs_genotype, sep='\t').iloc[:, 1:]
MHs_genotype_long = pd.melt(MHs_genotype_df, id_vars = 'IID', var_name = 'Locus.Allele', value_name = 'Genotype')
MHs_genotype_long[['Locus', 'Allele']] = MHs_genotype_long['Locus.Allele'].str.split('.', expand = True)
del(MHs_genotype_long['Locus.Allele'])
MHs_genotype_wide = MHs_genotype_long.pivot(index = ['IID', 'Allele'], columns = 'Locus', values = 'Genotype').reset_index().drop(columns = ['Allele'])

# Split MHs columns into SNPs columns
MHs_SNPs_genotype_df = MHs_genotype_wide.iloc[:, 0].copy()
for i in range(1, MHs_genotype_wide.shape[1]):
    if not IF_split_count_consistent(MHs_genotype_wide.iloc[:, i]):
        raise ValueError(f'SNP genotype is not consistent for column {MHs_genotype_wide.columns[i]}')
    MHs_SNPs_genotype_df = pd.concat([MHs_SNPs_genotype_df, Split_MHs_columns(MHs_genotype_wide.iloc[:, i])], axis = 1)

# Translate genotype to number
MHs_genotype_recode = MHs_SNPs_genotype_df.iloc[::2].copy().reset_index(drop = True)
for i in range(1, MHs_SNPs_genotype_df.shape[1]):
    major_allele = Get_major_allele(MHs_SNPs_genotype_df.iloc[:, i])
    MHs_genotype_recode.iloc[:, i] = Translate_genotype_to_number(MHs_SNPs_genotype_df.iloc[:, [0,i]], major_allele)

# Merge sample information to MHs_genotype_recode
sample_info_df = pd.read_csv(sample_info, sep='\t')
sample_info_df.columns = ['FID', 'IID', 'population']
MHs_genotype_recode_withInfo = pd.merge(sample_info_df[['IID', 'population']], MHs_genotype_recode, on = 'IID', how = 'right')

# Write MHs_SNPs_genotype_recode to file and run adegenet PCA
MHs_genotype_recode_withInfo.to_csv(output_prefix + '_MHs_SNPs_genotypes_for_adegenet_PCA.txt', index = False, sep = ' ')
_rscript = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rscript', 'PCA_using_adegenet.r')
cmd = f'Rscript {_rscript} -i {output_prefix}_MHs_SNPs_genotypes_for_adegenet_PCA.txt -o {output_prefix}'
os.system(cmd)
