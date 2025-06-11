# Description: This script contains functions that are loaded when the snakemake is initiated
# Modules required: numpy
# Date: 23 September 2024
# Author: Basile Pajot
#########################################################################################################################
########################  Functions   ###############################
######################## Get chromosome names  ###############################


def get_chromosome_name(path_chromosome_file):
    list_chromosomes = []
    with open(path_chromosome_file, "r") as f:
        for line in f.readlines():
            list_chromosomes.append(line.strip("\n"))
    return (list_chromosomes)
