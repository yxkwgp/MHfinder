import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = 'Stat allele frequency of MHs')
parser.add_argument('--mhs-genotype-file', '-g', type=str, required=True,
                       help='MHs genotype file')
parser.add_argument('--sample-info-file', '-s', type=str, required=True,
                       help='Sample information file')
parser.add_argument('--output-prefix', '-o', type=str, required=True,
                       help='Output filename prefix')
args = parser.parse_args()
mhs_genotype_file = args.mhs_genotype_file
sample_info_file = args.sample_info_file
output_prefix = args.output_prefix

#Reshape MH-genotype dataframe        
MHs_genotype_df = pd.read_csv(mhs_genotype_file, sep='\t').iloc[:, 1:]
MHs_genotype_long = pd.melt(MHs_genotype_df, id_vars = 'IID', var_name = 'Locus.Allele', value_name = 'Genotype')
MHs_genotype_long[['Locus', 'Allele']] = MHs_genotype_long['Locus.Allele'].str.split('.', expand = True)
del(MHs_genotype_long['Locus.Allele'])
MHs_genotype_wide = MHs_genotype_long.pivot(index = ['IID', 'Allele'], columns = 'Locus', values = 'Genotype').reset_index().drop(columns = ['Allele'])

# Add sample information
sample_info_df = pd.read_csv(sample_info_file, sep='\t').iloc[:, 1:]
sample_info_df.columns = ['IID', 'population']
MHs_genotype_wide_withInfo = pd.merge(MHs_genotype_wide, sample_info_df, on = 'IID', how = 'left')

# Stat allele frequency
dataPopOrder = MHs_genotype_wide_withInfo['population'].drop_duplicates().sort_values().tolist()
Locus_cols = [c for c in MHs_genotype_wide_withInfo.columns if c not in ['IID', 'population']]
MHs_genotype_long_for_freq = pd.melt(MHs_genotype_wide_withInfo, id_vars=['IID', 'population'], value_vars=Locus_cols, var_name='Locus', value_name='Allele')
MHs_genotype_long_for_freq = MHs_genotype_long_for_freq.dropna(subset=['Allele'])
allele_counts = MHs_genotype_long_for_freq.groupby(['Locus', 'population', 'Allele']).size().reset_index(name='n')
locus_pop_total = MHs_genotype_long_for_freq.groupby(['Locus', 'population']).size().reset_index(name='total')
allele_freq_long = allele_counts.merge(locus_pop_total, on=['Locus', 'population'], how='left')
allele_freq_long['freq'] = allele_freq_long['n'] / allele_freq_long['total']
allele_freq_table = allele_freq_long.pivot_table(index=['Locus', 'Allele'], columns='population', values='freq', fill_value=0).reset_index()
ordered_pop_cols = [pop for pop in dataPopOrder if pop in allele_freq_table.columns]
allele_freq_table = allele_freq_table[['Locus', 'Allele'] + ordered_pop_cols]
allele_freq_table.to_excel(output_prefix + '_alleleFreqTable.xlsx', index=False)
