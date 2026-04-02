import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = 'Filter MH-genotype file based on MHlist')
parser.add_argument('--mhs-genotype-file', '-g', type=str, required=True,
                       help='MHs genotype file')
parser.add_argument('--mhlist', '-m', type=str, required=True,
                       help='MHlist file')
parser.add_argument('--output', '-o', type=str, required=True,
                       help='Output filename')
args = parser.parse_args()
mhs_genotype_file = args.mhs_genotype_file
mhlist = args.mhlist
output = args.output

# Read in MHlist
MHslist_df = pd.read_csv(mhlist, header=None)
MHslist_df.columns = ['Locus']
MHslist_df['Locus.A1'] = MHslist_df['Locus'] + '.A1'
MHslist_df['Locus.A2'] = MHslist_df['Locus'] + '.A2'
MHs_colnames = MHslist_df[['Locus.A1', 'Locus.A2']].values.flatten().tolist()

# Read in MH-genotype file
MHs_genotype_df = pd.read_csv(mhs_genotype_file, sep='\t')
kept_colnames = MHs_genotype_df.columns.tolist()[:2] + MHs_colnames
MHs_genotype_df_filtered = MHs_genotype_df[kept_colnames]
MHs_genotype_df_filtered.to_csv(output, index=False, sep='\t')
