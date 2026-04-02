import argparse
import pandas as pd
import numpy as np
from sklearn import preprocessing

parser = argparse.ArgumentParser(description = 'Calculate Ae of MHs')
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

# Calculate Ae
MHs_genotype_wide_withInfo_group = MHs_genotype_wide_withInfo.groupby('population')
dataPopOrder = MHs_genotype_wide_withInfo['population'].drop_duplicates().tolist()
AeTable = np.zeros(shape = (MHs_genotype_wide_withInfo.shape[1] - 2, len(dataPopOrder)))
AeTable = pd.DataFrame(data = AeTable, index = MHs_genotype_wide.columns[1:], columns = dataPopOrder)
for pop in dataPopOrder:
    dataPart = MHs_genotype_wide_withInfo_group.get_group(pop).drop(['IID', 'population'], axis = 1)
    le = preprocessing.LabelEncoder()
    dataPartNum = [[1] * dataPart.shape[1]] * (dataPart.shape[0])
    dataPartNum = np.array(dataPartNum)
    for i in range(dataPart.shape[1]):
        le.fit(dataPart.iloc[:, i].drop_duplicates())
        dataPartNum[:, i] = le.transform(dataPart.iloc[:, i])
    for i in range(dataPartNum.shape[1]):
        alleleCount = np.bincount(dataPartNum[:,i])
        alleleFreq = alleleCount / dataPartNum.shape[0]
        Ae = 1 / np.sum(np.square(alleleFreq))
        AeTable.loc[MHs_genotype_wide.columns[1:][i], pop] = Ae

AeTable.index.name = None
AeTable.to_excel(output_prefix + '_AeTable.xlsx', index = True, header = True)
