# Libraries
import os
import pandas as pd
import numpy as np

########################  Functions   ###############################


def find_step_input(step, pop, raw_vcf_file, tmp_path):
    """
    This function is used to give different inputs for counting snps on the successfully filtered vcfs.

    Parameters:
    ------------------------------------
    step: str
        This is a short string (one of "Full_VCF", "Mac_data", "Removed_indels" or "Missing_data") to know at what 
        moment of the analysis we are. It will change with different vcf filtration steps.

    outputs_files: str
        This is the path to where the outputs are kept, including the temporary outputs that are used as input
        for several steps of the filtration process of vcfs

    final_output: str
        This is the path to where the permanent outputs are kept, including the Full VCF file that is used to do the filtration steps

    Returns:
    ------------------------------------
    Path (str): This is the path to use as input for the snp count.
    """
    if step == "1_Raw":
        return (raw_vcf_file)
    elif step == "2_Mac":
        return (tmp_path + pop + "/01_MAC/")
    elif step == "3_Biallelic":
        return (tmp_path + pop + "/02_Biall/")
    elif step == "4_Missing":
        return (tmp_path + pop + "/03_Missing/")
    elif step == "5_QUAL":
        return (tmp_path + pop + "/04_QUAL_filt/")
    elif step == "6_DP":
        return (tmp_path + pop + "/05_DP_filt/")
    elif step == "7_SP":
        return (tmp_path + pop + "/06_SP_filt/")
    elif step == "8_Maf":
        return (tmp_path + pop + "/07_Full_vcf/")
    elif step == "9_Hobs":
        return (tmp_path + pop + "/08_Hobs/")


def this_step(step):
    """
    This function is done to simplify the printing of the steps we are on for
    rules that are used several times.

    Parameters:
    ------------------------------------
    step: str
        This is a short string (one of "Full_VCF", "Mac_data", "Removed_indels" or "Missing_data") to know at what 
        moment of the analysis we are. It will change with different vcf filtration steps.

    Returns
    ------------------------------------
    : str
        Short string to print in the messages and the names of the rules
    """
    return (step.split('_')[1])
