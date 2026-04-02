import argparse
import os

parser = argparse.ArgumentParser(description = 'Maximum Likelihood Estimation of MHs genotypes using haplo.stats')
parser.add_argument('--plink-prefix', '-p', type=str, required=True,
                    help='PLINK binary file prefix')
parser.add_argument('--mhs-info-file', '-m', type=str, required=True,
                       help='MHs information file')
parser.add_argument('--output', '-o', type=str, required=True,
                       help='Output filename')
args = parser.parse_args()
plink_prefix = args.plink_prefix
mhs_info_file = args.mhs_info_file
output = args.output

# Get genotypes of SNPs in all candidate MHs

SNPs_in_MHs = set()
i = 1
infile = open(mhs_info_file, 'r')
for line in infile:
    if i % 3 == 0:
        line = line.strip().split('\t')
        SNPs_in_MHs.update(line)
    i += 1
infile.close()

outfile = open('SNPs_in_all_candidate_MHs.snplist', 'w')
for SNP in SNPs_in_MHs:
    outfile.write(SNP + '\n')
outfile.close()

cmd = f'plink --bfile {plink_prefix} --extract SNPs_in_all_candidate_MHs.snplist --recode --out SNPs_in_all_candidate_MHs'
os.system(cmd)

# Phase MHs using haplo.stats
_tools_dir = os.path.dirname(os.path.abspath(__file__))
_rscript = os.path.join(_tools_dir, 'rscript', 'MHs_phasing.R')
cmd = f'Rscript {_rscript} -p SNPs_in_all_candidate_MHs -m {mhs_info_file} -o {output}'
os.system(cmd)
