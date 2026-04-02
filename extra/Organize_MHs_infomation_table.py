import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description = 'Get basic information of MHs and organize them')
parser.add_argument('--mhs-info-file', '-m', type=str, required=True,
                       help='MHs information file')
parser.add_argument('--plink-prefix', '-p', type=str, required=True,
                       help='Plink binary file prefix')
parser.add_argument('--output-prefix', '-o', type=str, required=True,
                       help='Output filename prefix')
args = parser.parse_args()
mhs_info_file = args.mhs_info_file
plink_prefix = args.plink_prefix
output_prefix = args.output_prefix

# Get normal info and snplist
MHs_indexSNPs_dict = dict()
MHs_chr = dict()
MHs_unit = dict()
MHs_N_SNPs = dict()
MHs_Ae_list = []
SNPs = set()
infile = open(mhs_info_file, 'r')
i = 1
for line in infile:
    if i % 3 == 1:
        MH_name = line.strip().lstrip('>')
    elif i % 3 == 2:
        line = line.strip().split('\t')
        indexSNP = line[0].lstrip('*')
        MHs_indexSNPs_dict[MH_name] = indexSNP
        MHs_chr[MH_name] = int(line[1])
        MH_Ae = line[2].strip().split(', ')
        MH_Ae_df = pd.DataFrame([dict(x.split('=') for x in MH_Ae)])
        MH_Ae_df = MH_Ae_df.astype(float)
        MH_Ae_df['Name'] = MH_name
        MHs_Ae_list.append(MH_Ae_df)
    else:
        line = line.strip().split('\t')
        MHs_N_SNPs[MH_name] = len(line)
        SNPs.update(line)
        MHs_unit[MH_name] = '-'.join(line)
    i += 1
infile.close()
MHs_indexSNPs_df = pd.DataFrame(list(MHs_indexSNPs_dict.items()), columns = ['Name', 'Index SNP'])
MHs_chr_df = pd.DataFrame(list(MHs_chr.items()), columns = ['Name', 'Chr'])
MHs_unit_df = pd.DataFrame(list(MHs_unit.items()), columns = ['Name', 'Constitutional unit'])
MHs_N_SNPs_df = pd.DataFrame(list(MHs_N_SNPs.items()), columns = ['Name', 'N SNPs'])
MHs_Ae_df = pd.concat(MHs_Ae_list, ignore_index = True)

# Get bim file of all SNPs in MHs
snplist_file = f'{output_prefix}_all_SNPs.snplist'
outfile = open(snplist_file, 'w')
for SNP in SNPs:
    outfile.write(SNP + '\n')
outfile.close()
cmd = f'plink --bfile {plink_prefix} --extract {snplist_file} --make-bed --out {output_prefix}_all_SNPs'
os.system(cmd)

# Get BP
SNPs_bp = dict()
bim_file = open(f'{output_prefix}_all_SNPs.bim', 'r')
for line in bim_file:
    line = line.strip().split('\t')
    SNPs_bp[line[1]] = int(line[3])
bim_file.close()

# Get other info
MHs_MBP = dict()
MHs_length = dict()
infile = open(mhs_info_file, 'r')
i = 1
for line in infile:
    if i % 3 == 1:
        MH_name = line.strip().lstrip('>')
    elif i % 3 == 0:
        line = line.strip().split('\t')
        MHs_MBP[MH_name] = SNPs_bp[line[0]] / 1000000
        MHs_length[MH_name] = SNPs_bp[line[-1]] - SNPs_bp[line[0]] + 1
    i += 1
infile.close()
MHs_MBP_df = pd.DataFrame(list(MHs_MBP.items()), columns = ['Name', 'MBP'])
MHs_length_df = pd.DataFrame(list(MHs_length.items()), columns = ['Name', 'Length (bp)'])

# Merge all information and write to excel
MHs_organized_info_df = pd.merge(MHs_indexSNPs_df, MHs_chr_df, on = 'Name', how = 'left')
MHs_organized_info_df = pd.merge(MHs_organized_info_df, MHs_MBP_df, on = 'Name', how = 'left')
MHs_organized_info_df = pd.merge(MHs_organized_info_df, MHs_unit_df, on = 'Name', how = 'left')
MHs_organized_info_df = pd.merge(MHs_organized_info_df, MHs_N_SNPs_df, on = 'Name', how = 'left')
MHs_organized_info_df = pd.merge(MHs_organized_info_df, MHs_length_df, on = 'Name', how = 'left')
MHs_organized_info_df = pd.merge(MHs_organized_info_df, MHs_Ae_df, on = 'Name', how = 'left')
MHs_organized_info_df.to_excel(f'{output_prefix}_organized_MHs_info.xlsx', index = False)
