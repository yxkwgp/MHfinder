import argparse
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.naive_bayes import CategoricalNB
from scipy.special import logsumexp
from sklearn import metrics
import os

parser = argparse.ArgumentParser(description = 'Leave-one-out cross-validation of MHs')
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

def joint_log_likelihood(X):
    jll = np.zeros(CNB.class_count_.shape[0])
    for j in range(CNB.n_features_in_):
        indice1 = X[0, j]
        indice2 = X[1, j]
        if indice1 < CNB.feature_log_prob_[j].shape[1]:
            jll += CNB.feature_log_prob_[j][:, indice1]
        if indice2 < CNB.feature_log_prob_[j].shape[1]:
            jll += CNB.feature_log_prob_[j][:, indice2]
    total_ll = jll + CNB.class_log_prior_
    return total_ll

# Reshape MH-genotype dataframe        
MHs_genotype_df = pd.read_csv(mhs_genotype_file, sep='\t').iloc[:, 1:]
MHs_genotype_long = pd.melt(MHs_genotype_df, id_vars = 'IID', var_name = 'Locus.Allele', value_name = 'Genotype')
MHs_genotype_long[['Locus', 'Allele']] = MHs_genotype_long['Locus.Allele'].str.split('.', expand = True)
del(MHs_genotype_long['Locus.Allele'])
MHs_genotype_wide = MHs_genotype_long.pivot(index = ['IID', 'Allele'], columns = 'Locus', values = 'Genotype').reset_index().drop(columns = ['Allele'])

# Add sample information
sample_info_df = pd.read_csv(sample_info_file, sep='\t').iloc[:, 1:]
sample_info_df.columns = ['IID', 'population']
MHs_genotype_wide_withInfo = pd.merge(MHs_genotype_wide, sample_info_df, on = 'IID', how = 'left')

# Prepare y for LOOCV
populationTrue = list(MHs_genotype_wide_withInfo['population'])[::2]
popEncoder = preprocessing.LabelEncoder()
popEncoder.fit(MHs_genotype_wide_withInfo['population'].drop_duplicates())
yTwice = popEncoder.transform(MHs_genotype_wide_withInfo['population'])

# Label encoding MH-genotype
le = preprocessing.LabelEncoder()
dataNum = [[1] * (MHs_genotype_wide_withInfo.shape[1] - 2)] * MHs_genotype_wide_withInfo.shape[0]
dataNum = np.array(dataNum)
for i in range(1, MHs_genotype_wide_withInfo.shape[1] - 1):
    le.fit(MHs_genotype_wide_withInfo.iloc[:, i].drop_duplicates())
    dataNum[:, i - 1] = le.transform(MHs_genotype_wide_withInfo.iloc[:, i])

# Build model and predict
probTable = np.zeros(shape=(len(populationTrue), popEncoder.classes_.shape[0]))
for i in range(len(populationTrue)):
    trainX = np.delete(dataNum, [2 * i, 2 * i + 1], axis=0)
    testX = dataNum[2 * i: 2 * i + 2,]
    trainY = np.concatenate([yTwice[:2 * i], yTwice[2 * i + 2:]])
    trainY = np.array(trainY)
    CNB = CategoricalNB()
    CNB.fit(trainX, trainY)
    jll = joint_log_likelihood(testX)
    log_prob_all = logsumexp(jll)
    normalizedLogProb = jll - log_prob_all
    probTable[i,] = np.exp(normalizedLogProb)
yPred = CNB.classes_[np.argmax(probTable, axis=1)]
populationPred = popEncoder.inverse_transform(yPred)

# Calculate evaluation metrics
## Accuracy
Accuracy = metrics.accuracy_score(list(populationTrue), list(populationPred))
accFile = open(f'{output_prefix}_accuracy.txt', 'w')
accFile.write(f'{Accuracy}\n')
accFile.close()
## Confusion matrix
Confusion_matrix = metrics.confusion_matrix(list(populationTrue), list(populationPred)).T
Confusion_matrix_labels = popEncoder.inverse_transform(np.arange(Confusion_matrix.shape[0]))
Confusion_matrix_df = pd.DataFrame(Confusion_matrix, index = Confusion_matrix_labels, columns = Confusion_matrix_labels)
Confusion_matrix_df.index.name = 'Prediction'
Confusion_matrix_df.columns = pd.MultiIndex.from_arrays([
    ['Target'] * len(Confusion_matrix_labels),
    Confusion_matrix_labels
])
Confusion_matrix_df.to_excel(f'{output_prefix}_confusion_matrix.xlsx')
## AUC, Specificity, Sensitivity, PPV and NPV
lb = preprocessing.LabelBinarizer()
populationTrueforROC = lb.fit_transform(populationTrue)
poulationPredforROC = lb.fit_transform(populationPred)
AUCs = np.zeros(popEncoder.classes_.shape[0])
Specificities = np.zeros(popEncoder.classes_.shape[0])
Sensitivities = np.zeros(popEncoder.classes_.shape[0])
PPVs = np.zeros(popEncoder.classes_.shape[0])
NPVs = np.zeros(popEncoder.classes_.shape[0])
for i in range(popEncoder.classes_.shape[0]):
    AUCs[i] = metrics.roc_auc_score(populationTrueforROC[:,i], probTable[:,i])
    Sensitivities[i] = metrics.recall_score(populationTrueforROC[:,i], poulationPredforROC[:,i], pos_label = 1)
    Specificities[i] = metrics.recall_score(populationTrueforROC[:,i], poulationPredforROC[:,i], pos_label = 0)
    PPVs[i] = metrics.precision_score(populationTrueforROC[:,i], poulationPredforROC[:,i], pos_label = 1)
    NPVs[i] = metrics.precision_score(populationTrueforROC[:,i], poulationPredforROC[:,i], pos_label = 0)
evaluationMetrics = pd.DataFrame(np.column_stack([AUCs, Specificities, Sensitivities, PPVs, NPVs]), columns = ['AUC', 'Specificity', 'Sensitivity', 'PPV', 'NPV'])
evaluationMetrics.index = Confusion_matrix_labels
evaluationMetrics.to_excel(f'{output_prefix}_evaluation_metrics.xlsx')

# Draw confusion matrix plot
## Write yPred and yTrue to file
yPredFile = open(f'{output_prefix}_yPred.txt', 'w')
yTrueFile = open(f'{output_prefix}_yTrue.txt', 'w')
for i in range(len(populationPred)):
    yPredFile.write(f'{populationPred[i]}\n')
    yTrueFile.write(f'{populationTrue[i]}\n')
yPredFile.close()
yTrueFile.close()
## Draw confusion matrix plot
_rscript = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rscript', 'Draw_confusion_matrix.r')
cmd = f'Rscript {_rscript} -p {output_prefix}_yPred.txt -t {output_prefix}_yTrue.txt -o {output_prefix}_confusion_matrix.tiff'
os.system(cmd)
