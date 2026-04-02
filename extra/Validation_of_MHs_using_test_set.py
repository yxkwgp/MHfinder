import argparse
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.naive_bayes import CategoricalNB
from scipy.special import logsumexp
from sklearn import metrics
import os

parser = argparse.ArgumentParser(description = 'Validation of MHs using test set')
parser.add_argument('--training-genotype-file', '-g', type=str, required=True,
                       help='Training set genotype file')
parser.add_argument('--test-genotype-file', '-t', type=str, required=True,
                       help='Test set genotype file')
parser.add_argument('--sample-info-file', '-s', type=str, required=True,
                       help='Sample information file')
parser.add_argument('--output-prefix', '-o', type=str, required=True,
                       help='Output filename prefix')
args = parser.parse_args()
training_genotype_file = args.training_genotype_file
test_genotype_file = args.test_genotype_file
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

# Read in sample information
sample_info_df = pd.read_csv(sample_info_file, sep='\t').iloc[:, 1:]
sample_info_df.columns = ['IID', 'population']

# Read in and reshape training set genotype     
training_genotype_df = pd.read_csv(training_genotype_file, sep='\t').iloc[:, 1:]
training_genotype_long = pd.melt(training_genotype_df, id_vars = 'IID', var_name = 'Locus.Allele', value_name = 'Genotype')
training_genotype_long[['Locus', 'Allele']] = training_genotype_long['Locus.Allele'].str.split('.', expand = True)
del(training_genotype_long['Locus.Allele'])
training_genotype_wide = training_genotype_long.pivot(index = ['IID', 'Allele'], columns = 'Locus', values = 'Genotype').reset_index().drop(columns = ['Allele'])
training_genotype_wide_withInfo = pd.merge(training_genotype_wide, sample_info_df, on = 'IID', how = 'left')
populationTrueTrain = list(training_genotype_wide_withInfo['population'])[::2]

# Read in and reshape test set genotype
test_genotype_df = pd.read_csv(test_genotype_file, sep='\t').iloc[:, 1:]
test_genotype_long = pd.melt(test_genotype_df, id_vars = 'IID', var_name = 'Locus.Allele', value_name = 'Genotype')
test_genotype_long[['Locus', 'Allele']] = test_genotype_long['Locus.Allele'].str.split('.', expand = True)
del(test_genotype_long['Locus.Allele'])
test_genotype_wide = test_genotype_long.pivot(index = ['IID', 'Allele'], columns = 'Locus', values = 'Genotype').reset_index().drop(columns = ['Allele'])
test_genotype_wide_withInfo = pd.merge(test_genotype_wide, sample_info_df, on = 'IID', how = 'left')
populationTrueTest = list(test_genotype_wide_withInfo['population'])[::2]

# Merge training and test set genotype, add sample information
dataMerge_withInfo = pd.concat([training_genotype_wide_withInfo, test_genotype_wide_withInfo], ignore_index = True)

# Prepare y
popEncoder = preprocessing.LabelEncoder()
popEncoder.fit(dataMerge_withInfo['population'].drop_duplicates())
yTrainTwice = popEncoder.transform(training_genotype_wide_withInfo['population'])

# LabelEncoder coding
le = preprocessing.LabelEncoder()
dataNum = [[1] * (dataMerge_withInfo.shape[1] - 2)] * dataMerge_withInfo.shape[0]
dataNum = np.array(dataNum)
for i in range(1, dataMerge_withInfo.shape[1] - 1):
    le.fit(dataMerge_withInfo.iloc[:, i].drop_duplicates())
    dataNum[:, i - 1] = le.transform(dataMerge_withInfo.iloc[:, i])

# Build model
trainX = dataNum[:training_genotype_wide_withInfo.shape[0],]
trainY = np.array(yTrainTwice)
CNB = CategoricalNB()
CNB.fit(trainX, trainY)

# Predict test set
probTable = np.zeros(shape=(len(populationTrueTest), popEncoder.classes_.shape[0]))
for i in range(len(populationTrueTest)):
    testX = dataNum[2 * i + training_genotype_wide_withInfo.shape[0]: 2 * i + 2 + training_genotype_wide_withInfo.shape[0],]
    jll = joint_log_likelihood(testX)
    log_prob_all = logsumexp(jll)
    normalizedLogProb = jll - log_prob_all
    probTable[i,] = np.exp(normalizedLogProb)
yPred = CNB.classes_[np.argmax(probTable, axis=1)]
yPred = list(yPred)
populationPredTest = popEncoder.inverse_transform(yPred)

# Calculate evaluation metrics
## Accuracy
Accuracy = metrics.accuracy_score(list(populationTrueTest), list(populationPredTest))
accFile = open(f'{output_prefix}_accuracy.txt', 'w')
accFile.write(f'{Accuracy}\n')
accFile.close()
## Confusion matrix
Confusion_matrix = metrics.confusion_matrix(list(populationTrueTest), list(populationPredTest)).T
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
populationTrueforROC = lb.fit_transform(populationTrueTest)
poulationPredforROC = lb.fit_transform(populationPredTest)
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
for i in range(len(populationPredTest)):
    yPredFile.write(f'{populationPredTest[i]}\n')
    yTrueFile.write(f'{populationTrueTest[i]}\n')
yPredFile.close()
yTrueFile.close()
## Draw confusion matrix plot
_rscript = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rscript', 'Draw_confusion_matrix.r')
cmd = f'Rscript {_rscript} -p {output_prefix}_yPred.txt -t {output_prefix}_yTrue.txt -o {output_prefix}_confusion_matrix.tiff'
os.system(cmd)
