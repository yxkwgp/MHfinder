import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = 'Prepare MH-genotype file for Infocalc')
parser.add_argument('--mhs-genotype-file', '-g', type=str, required=True,
                       help='MHs genotype file')
parser.add_argument('--sample-info', '-s', type=str, required=True,
                    help='Sample information file')
parser.add_argument('--output', '-o', type=str, required=True,
                       help='Output filename')
args = parser.parse_args()
mhs_genotype_file = args.mhs_genotype_file
sample_info = args.sample_info
output = args.output

def check_split_count_consistent(pd_series):
    split_counts = pd_series.str.split('-').str.len()
    unique_counts = split_counts.unique()
    if len(unique_counts) == 1:
        return True
    else:
        return False

def translate_genotype_to_number(pd_series):
    pd_series = pd_series.str.split('-')
    SNP_num = len(pd_series[0])
    translated_genotype = pd.DataFrame(index=pd_series.index, columns=range(SNP_num))
    for i in range(SNP_num):
        SNP_order = pd_series.str[i].value_counts().index.tolist()
        SNP_order_dict = {value: str(index + 1) for index, value in enumerate(SNP_order)}
        translated_genotype.iloc[:, i] = pd_series.str[i].map(SNP_order_dict)
    translated_genotype_result = translated_genotype.agg(''.join, axis=1)
    return translated_genotype_result

# Reshape MH-genotype dataframe        
MHs_genotype_df = pd.read_csv(mhs_genotype_file, sep='\t').iloc[:, 1:]
MHs_genotype_long = pd.melt(MHs_genotype_df, id_vars = 'IID', var_name = 'Locus.Allele', value_name = 'Genotype')
MHs_genotype_long[['Locus', 'Allele']] = MHs_genotype_long['Locus.Allele'].str.split('.', expand = True)
del(MHs_genotype_long['Locus.Allele'])
MHs_genotype_wide = MHs_genotype_long.pivot(index = ['IID', 'Allele'], columns = 'Locus', values = 'Genotype').reset_index().drop(columns = ['Allele'])

# Add sample information
sample_info_df = pd.read_csv(sample_info, sep='\t')
sample_info_df.columns = ['FID', 'IID', 'population']
MHs_genotype_wide_withInfo = pd.merge(MHs_genotype_wide, sample_info_df, on = 'IID', how = 'left')
MHs_genotype_wide_withInfo['pop1'] = '-9'
MHs_genotype_wide_withInfo['pop2'] = '-9'
front_cols = ['IID', 'FID', 'population', 'pop1', 'pop2']
other_cols = [colname for colname in MHs_genotype_wide_withInfo.columns if colname not in front_cols]
MHs_genotype_wide_withInfo = MHs_genotype_wide_withInfo[front_cols + other_cols]

# Translate genotype to number
for i in range(MHs_genotype_wide_withInfo.shape[1] - 5):
    if not check_split_count_consistent(MHs_genotype_wide_withInfo.iloc[:, i + 5]):
        raise ValueError(f'SNP genotype is not consistent for column {MHs_genotype_wide_withInfo.columns[i + 5]}')
    MHs_genotype_wide_withInfo.iloc[:, i + 5] = translate_genotype_to_number(MHs_genotype_wide_withInfo.iloc[:, i + 5])

# Write output file
output_header = MHs_genotype_wide_withInfo.columns[5:].tolist()
output_rows = MHs_genotype_wide_withInfo.values.tolist()
outfile = open(output, 'w')
outfile.write('\t'.join(str(c) for c in output_header) + '\n')
for row in output_rows:
    outfile.write('\t'.join(str(x) for x in row) + '\n')
outfile.close()
