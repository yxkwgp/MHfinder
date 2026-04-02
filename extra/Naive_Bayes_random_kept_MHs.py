import argparse
from math import comb
from random import seed
from random import sample
from itertools import combinations
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.naive_bayes import CategoricalNB
from scipy.special import logsumexp
from sklearn import metrics
from os import makedirs
from os import system
import multiprocessing as mp
import time

parser = argparse.ArgumentParser(description = 'Evaluation of random kept MHs using Naive Bayes')
parser.add_argument('--mhs-genotype-file', '-g', type=str, required=True,
                       help='MHs genotype file')
parser.add_argument('--sample-info-file', '-s', type=str, required=True,
                       help='Sample information file')
parser.add_argument('--random-seed', '-r', type = int, default = 123,
                       help='Random seed')
parser.add_argument('--times', '-t', type=int, default=10000,
                       help='Number of random subset simulations for each kept size (default: 10000)')
parser.add_argument('--output-dir', '-d', type=str, required=True,
                       help='Output directory')
args = parser.parse_args()
if args.times < 1:
    parser.error('--times must be >= 1')
mhs_genotype_file = args.mhs_genotype_file
sample_info_file = args.sample_info_file
random_seed = args.random_seed
sim_times = args.times
output_dir = args.output_dir
output_dir = output_dir.rstrip('/') + '/'

def generate_random_combinations(total_elements, sample_size, num_combinations, random_seed=None):
    if random_seed is not None:
        seed(random_seed)
        np.random.seed(random_seed)
    
    num_all_combinations = comb(total_elements, sample_size)
    if num_combinations > num_all_combinations:
        random_combinations = list(combinations(list(range(total_elements)), sample_size))
        random_combinations = [set(element) for element in random_combinations]
    elif num_combinations > num_all_combinations * 0.05:
        all_combinations = list(combinations(list(range(total_elements)), sample_size))
        all_combinations = [set(element) for element in all_combinations]
        random_combinations = sample(all_combinations, num_combinations)
    else:
        random_combinations = []
        while len(random_combinations) < num_combinations:
            elements = list(range(total_elements))
            combination = set(sample(elements, sample_size))
            if combination not in random_combinations:
                random_combinations.append(combination)
    return random_combinations

def joint_log_likelihood(X, CNB):
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

def single_simulation(simulation_data):
    # Unpack data
    j = simulation_data['simulation_id']
    random_combination = simulation_data['random_combination']
    dataNum = simulation_data['dataNum']
    populationTrue = simulation_data['populationTrue']
    yTwice = simulation_data['yTwice']
    popEncoder = simulation_data['popEncoder']
    populationTrueforROC = simulation_data['populationTrueforROC']
    n_MHs = simulation_data['n_MHs']
    # Run single simulation
    try:
        dataNum_degraded = dataNum[:, list(random_combination)]
        MHs_kept_IF = np.isin(list(range(n_MHs)), list(random_combination)).astype(int)
        probTable = np.zeros(shape=(len(populationTrue), popEncoder.classes_.shape[0]))
        
        for k in range(len(populationTrue)):
            trainX = np.delete(dataNum_degraded, [2 * k, 2 * k + 1], axis=0)
            testX = dataNum_degraded[2 * k: 2 * k + 2,]
            trainY = np.concatenate([yTwice[:2 * k], yTwice[2 * k + 2:]])
            trainY = np.array(trainY)
            CNB = CategoricalNB()
            CNB.fit(trainX, trainY)
            jll = joint_log_likelihood(testX, CNB)
            log_prob_all = logsumexp(jll)
            normalizedLogProb = jll - log_prob_all
            probTable[k,] = np.exp(normalizedLogProb)
        
        yPred = CNB.classes_[np.argmax(probTable, axis=1)]
        populationPred = list(popEncoder.inverse_transform(yPred))
        accuracy = metrics.accuracy_score(populationTrue, populationPred)
        
        AUC = []
        for l in range(popEncoder.classes_.shape[0]):
            AUC.append(metrics.roc_auc_score(populationTrueforROC[:,l], probTable[:,l]))
        AUC = np.array(AUC)
        AUC_mean = AUC.mean()
        AUC_all = np.insert(AUC, 0, AUC_mean)
        
        return {
            'simulation_id': j,
            'accuracy': accuracy,
            'AUC_all': AUC_all,
            'populationPred': populationPred,
            'MHs_kept_IF': MHs_kept_IF,
            'success': True
        }
        
    except Exception as e:
        print(f'Simulation {j} failed: {e}')
        return {
            'simulation_id': j,
            'accuracy': 0.0,
            'AUC_all': np.zeros(popEncoder.classes_.shape[0] + 1),
            'populationPred': [None] * len(populationTrue),
            'MHs_kept_IF': np.zeros(n_MHs),
            'success': False,
            'error': str(e)
        }

def run_parallel_simulations(simulation_data_list, num_processes=10):
    start_time = time.time()
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(single_simulation, simulation_data_list)
    end_time = time.time()
    print(f'Parallel simulations completed, time: {end_time - start_time:.2f} seconds')
    # Sort by simulation_id to ensure consistency
    results.sort(key=lambda x: x['simulation_id'])
    return results

if __name__ == '__main__':
    # Set random seed for reproducibility
    RANDOM_SEED = random_seed
    
    # Set multiprocessing start method
    mp.set_start_method('spawn', force=True)
    
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

    # Prepare y
    populationTrue = list(MHs_genotype_wide_withInfo['population'])[::2]
    popEncoder = preprocessing.LabelEncoder()
    popEncoder.fit(MHs_genotype_wide_withInfo['population'].drop_duplicates())
    yTwice = popEncoder.transform(MHs_genotype_wide_withInfo['population'])
    lb = preprocessing.LabelBinarizer()
    populationTrueforROC = lb.fit_transform(populationTrue)
    IID = list(MHs_genotype_wide_withInfo['IID'])[::2]

    # Label encoding MH-genotype
    le = preprocessing.LabelEncoder()
    dataNum = [[1] * (MHs_genotype_wide_withInfo.shape[1] - 2)] * MHs_genotype_wide_withInfo.shape[0]
    dataNum = np.array(dataNum)
    for i in range(1, MHs_genotype_wide_withInfo.shape[1] - 1):
        le.fit(MHs_genotype_wide_withInfo.iloc[:, i].drop_duplicates())
        dataNum[:, i - 1] = le.transform(MHs_genotype_wide_withInfo.iloc[:, i])

    # Set number of processes based on max number of processes
    num_processes = min(10, mp.cpu_count() - 1)
    print(f"Use {num_processes} processes for parallel computation")

    # Simulate degradation
    n_MHs = dataNum.shape[1]
    print(f"Start simulating degradation, total features: {n_MHs}")

    for i in range(n_MHs - 1, 0, -1):
        print(f"\n=== Processing {i} features ===")
        
        # Generate random combinations
        print("Generate random combinations...")
        random_combination = generate_random_combinations(n_MHs, i, sim_times, random_seed=RANDOM_SEED)
        
        # Create output directory
        makedirs(f'{output_dir}/kept_{i}', exist_ok = True)
        
        # Prepare simulation data
        print("Prepare simulation data...")
        simulation_data_list = []
        for j in range(len(random_combination)):
            simulation_data = {
                'simulation_id': j,
                'random_combination': random_combination[j],
                'dataNum': dataNum,
                'populationTrue': populationTrue,
                'yTwice': yTwice,
                'popEncoder': popEncoder,
                'populationTrueforROC': populationTrueforROC,
                'n_MHs': n_MHs
            }
            simulation_data_list.append(simulation_data)
        
        # Run parallel simulations
        print(f"Start {len(random_combination)} parallel simulations...")
        results = run_parallel_simulations(simulation_data_list, num_processes)
        
        # Collect and save results
        print("Collect and save results...")
        ACC_file = open(f'{output_dir}/kept_{i}/ACC.txt', 'w')
        AUC_file = open(f'{output_dir}/kept_{i}/AUC.txt', 'w')
        populationPredAll = []
        MHs_kept_IF_all = []
        
        successful_simulations = 0
        for result in results:
            if result['success']:
                # Write to files
                ACC_file.write(str(result['accuracy']) + '\n')
                AUC_all_str = '\t'.join(map(str, list(result['AUC_all'])))
                AUC_file.write(AUC_all_str + '\n')
                
                # Collect data
                populationPredAll.append(result['populationPred'])
                MHs_kept_IF_all.append(result['MHs_kept_IF'])
                successful_simulations += 1
            else:
                print(f"Simulation {result['simulation_id']} failed: {result.get('error', 'Unknown error')}")
        
        ACC_file.close()
        AUC_file.close()
        
        print(f'{successful_simulations} times simulations are successful')
        
        # Save DataFrames
        if populationPredAll:
            populationPredAll_df = pd.DataFrame(populationPredAll, columns = IID)
            populationPredAll_df.to_csv(f'{output_dir}/kept_{i}/populationPredAll.txt', sep = '\t', index = False)
            
            MHs_kept_IF_all_df = pd.DataFrame(MHs_kept_IF_all, columns = MHs_genotype_wide.columns[1:])
            MHs_kept_IF_all_df.to_csv(f'{output_dir}/kept_{i}/MHs_kept_IF_all.txt', sep = '\t', index = False)
            
            print(f'Results saved to {output_dir}/kept_{i}/ directory')
        else:
            print(f'Warning: no successful simulation results')
