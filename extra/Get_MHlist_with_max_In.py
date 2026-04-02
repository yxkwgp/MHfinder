import argparse
import pandas as pd

parser = argparse.ArgumentParser(description = 'Reduce loci with the same key SNP to one based on In')
parser.add_argument('--In-file', '-i', type=str, required=True,
                       help='In file')
parser.add_argument('--output', '-o', type=str, required=True,
                       help='Output filename')
args = parser.parse_args()
In_file = args.In_file
output = args.output

# Read in In file
df = pd.read_csv(In_file, sep=r'\s+', skipfooter=2, engine='python')
df['population'] = df['Locus'].str.split('_').str[0]
df['key_SNP_index'] = df['Locus'].str.split('_').str[1].astype(int)
df['MH_index'] = df['Locus'].str.split('_').str[2].astype(int)
index_maxIn = df.groupby(['population', 'key_SNP_index'])['I_n'].idxmax()
df_maxIn = df.loc[index_maxIn]
df_maxIn_sorted = df_maxIn.sort_values(by=['population', 'key_SNP_index']).reset_index(drop=True)
df_maxIn_sorted['Locus'].to_csv(output, index = False, header = False)
