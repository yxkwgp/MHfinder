import os
import subprocess
import sys


def _run_tool(script_name):
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), script_name)
    cmd = [sys.executable, script_path] + sys.argv[1:]
    completed = subprocess.run(cmd)
    raise SystemExit(completed.returncode)


def mle_mh_genotypes_cli():
    _run_tool("MLE_for_MHs_genotypes.py")


def plot_scgwas_manhattan_qq_cli():
    _run_tool("Draw_manhattan_qq_of_SCGWAS.py")


def plot_mh_pca_cli():
    _run_tool("Draw_PCA_of_MHs.py")


def loocv_mhs_cli():
    _run_tool("LOOCV_of_MHs.py")


def validate_mhs_testset_cli():
    _run_tool("Validation_of_MHs_using_test_set.py")


def organize_mh_info_table_cli():
    _run_tool("Organize_MHs_infomation_table.py")


def stat_mh_allele_freq_cli():
    _run_tool("Stat_allele_freq_of_MHs.py")


def stat_mh_indexsnps_cli():
    _run_tool("Stat_MHs_indexSNPs_info.py")


def calc_mh_ae_cli():
    _run_tool("Calculate_Ae_of_MHs.py")


def naivebayes_random_kept_mhs_cli():
    _run_tool("Naive_Bayes_random_kept_MHs.py")


def filter_mh_genotype_by_list_cli():
    _run_tool("Filter_MH_genotype_file_based_MHlist.py")


def get_mhlist_max_in_cli():
    _run_tool("Get_MHlist_with_max_In.py")


def prep_infocalc_input_cli():
    _run_tool("Prep_file_for_Infocalc.py")


def auc_acc_decline_with_n_mhs_cli():
    _tools_dir = os.path.dirname(os.path.abspath(__file__))
    r_script = os.path.join(_tools_dir, "AUC_ACC_decline_with_n_MHs.r")
    cmd = ["Rscript", r_script] + sys.argv[1:]
    completed = subprocess.run(cmd)
    raise SystemExit(completed.returncode)

