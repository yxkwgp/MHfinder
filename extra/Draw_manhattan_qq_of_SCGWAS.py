import argparse
import os

parser = argparse.ArgumentParser(description = 'Draw Manhattan and QQ plots of SCGWAS')
parser.add_argument('--mhs-info-file', '-m', type=str, required=True,
                       help='MHs information file')
parser.add_argument('--p-value-threshold', '-p', type=float, required=True,
                       help='P-value threshold of SCGWAS')
parser.add_argument('--scgwas-dir', '-d', type=str, required=True,
                       help='directory of SCGWAS')
parser.add_argument('--output-prefix', '-o', type=str, required=True,
                       help='Output filename prefix')
args = parser.parse_args()
mhs_info_file = args.mhs_info_file
p_value_threshold = args.p_value_threshold
scgwas_dir = args.scgwas_dir
output_prefix = args.output_prefix

scgwas_dir = scgwas_dir.rstrip('/') + '/'

# Prepare .assoc.logistic files during reading MHs information file
infile = open(mhs_info_file, 'r')

## Initialize the first population
line1 = infile.readline().strip().lstrip('>').split('_')
popAll = []
pop = line1[0]
popAll.append(pop)
roundTime = int(line1[1])
roundTimeRem = roundTime
output_dir = os.path.dirname(output_prefix)
if output_dir != '':
    output_dir = output_dir.rstrip('/') + '/'
outfile_name = f'{output_dir}SCGWAS_{pop}_merge.assoc.logistic'
outfile = open(outfile_name, 'w')
outfile.write(' CHR              SNP         BP   A1       TEST    NMISS         OR       SE      L95      U95         STAT            P \n')
line2 = infile.readline().strip().split('\t')
indexSNP = line2[0].lstrip('*')
assoc_file_name = f'{scgwas_dir}SCGWAS_{pop}_{roundTime}.assoc.logistic'
assoc_file = open(assoc_file_name, 'r')
for gwas_line in assoc_file:
    gwas_line_tmp = gwas_line.strip().split()
    if gwas_line_tmp[1] == indexSNP:
        outfile.write(gwas_line)
        break
assoc_file.close()
infile.readline()
i = 1
for line in infile:

    ## Read in a MH and its corresponding pop and round
    if i % 3 == 1:
        line = line.strip().lstrip('>').split('_')
        pop = line[0]
        roundTime = int(line[1])
    elif i % 3 == 2:
        line = line.strip().split('\t')
        indexSNP = line[0].lstrip('*')

        ## Get SNPs need to be deleted from the background gwas results
        if pop != popAll[-1]:
            final_gwas_filtered_file = f'{scgwas_dir}SCGWAS_{popAll[-1]}_{roundTimeRem + 1}.filter.assoc.logistic'
            final_gwas_filtered_file = open(final_gwas_filtered_file, 'r')
            del_SNPs = []
            for gwas_line in final_gwas_filtered_file:
                gwas_line = gwas_line.strip().split()
                del_SNPs.append(gwas_line[0])
            final_gwas_filtered_file.close()

            ## Get background SNPs from the final round gwas of this pop
            final_gwas_unfiltered_file = f'{scgwas_dir}SCGWAS_{popAll[-1]}_{roundTimeRem + 1}.assoc.logistic'
            final_gwas_unfiltered_file = open(final_gwas_unfiltered_file, 'r')
            final_gwas_unfiltered_file.readline()
            for gwas_line in final_gwas_unfiltered_file:
                gwas_line_tmp = gwas_line.strip().split()
                if gwas_line_tmp[1] not in del_SNPs:
                    outfile.write(gwas_line)
            final_gwas_unfiltered_file.close()
            outfile.close()

            ## Initialize the next population
            outfile_name = f'{output_dir}SCGWAS_{pop}_merge.assoc.logistic'
            outfile = open(outfile_name, 'w')
            outfile.write(' CHR              SNP         BP   A1       TEST    NMISS         OR       SE      L95      U95         STAT            P \n')
            popAll.append(pop)
        assoc_file_name = f'{scgwas_dir}SCGWAS_{pop}_{roundTime}.assoc.logistic'
        assoc_file = open(assoc_file_name, 'r')
        for gwas_line in assoc_file:
            gwas_line_tmp = gwas_line.strip().split()
            if gwas_line_tmp[1] == indexSNP:
                outfile.write(gwas_line)
                break
        assoc_file.close()
        roundTimeRem = roundTime
    i += 1

## Background of final population
final_gwas_filtered_file = f'{scgwas_dir}SCGWAS_{popAll[-1]}_{roundTimeRem + 1}.filter.assoc.logistic'
final_gwas_filtered_file = open(final_gwas_filtered_file, 'r')
del_SNPs = []
for gwas_line in final_gwas_filtered_file:
    gwas_line = gwas_line.strip().split()
    del_SNPs.append(gwas_line[0])
final_gwas_filtered_file.close()
final_gwas_unfiltered_file = f'{scgwas_dir}SCGWAS_{popAll[-1]}_{roundTimeRem + 1}.assoc.logistic'
final_gwas_unfiltered_file = open(final_gwas_unfiltered_file, 'r')
final_gwas_unfiltered_file.readline()
for gwas_line in final_gwas_unfiltered_file:
    gwas_line_tmp = gwas_line.strip().split()
    if gwas_line_tmp[1] not in del_SNPs:
        outfile.write(gwas_line)
final_gwas_unfiltered_file.close()
outfile.close()

# End of MHs information and run Rscript
infile.close()
gwas_result_files = [f'{output_dir}SCGWAS_{pop}_merge.assoc.logistic' for pop in popAll]
_rscript = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rscript', 'Manhattan_qq_plot_of_many_gwas_result.r')
gwas_g = ','.join(gwas_result_files)
cmd = f'Rscript {_rscript} -p {p_value_threshold} -g {gwas_g} -o {output_prefix}'
os.system(cmd)
