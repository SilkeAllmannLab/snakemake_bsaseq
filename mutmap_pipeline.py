#!/usr/bin/env python3

import yaml
import pandas as pd
import os
import subprocess

##################
# Helper functions
##################
def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples_df.columns:
        return True
    else:
        return pd.isnull(samples_df.loc[(sample), "fq2"])

def get_trimmed_files_per_sample_type(trimmed_dir = "trimmed", sample_type = "cultivar"):
    """
    This function:
      1. Collects the files names of the trimmed files inside the trimmed directory
      2. Filter to keep only the "cultivar" or "bulk" trimmed files.
      3. Returns a list with one or two elements. 
    """
    if os.path.isdir(trimmed_dir):
        pass
    else:
        print("Directory with trimmed file does not exist")

    trimmed_files = os.listdir(trimmed_dir) # collects all trimmed files
    trimmed_files.sort()                    # to have R1 before R2

    cultivar_files = [trimmed_dir + f for f in trimmed_files if sample_type in f]
    bulk_files = [trimmed_dir + f for f in trimmed_files if sample_type in f]
    if sample_type == "cultivar":
        return ",".join(cultivar_files) # fastq of cultivar. If you specify fastq, please separate pairs by comma, e.g. -b fastq1,fastq2.
    elif sample_type == "bulk":
        return ",".join(bulk_files)     # fastq of mutant bulk. If you specify fastq, please separate pairs by comma, e.g. -b fastq1,fastq2.
    else:
        "Please check that sample_type == 'cultivar' or sample_type == 'bulk'"

############################
# Import pipeline parameters
############################
with open("config/config.yaml") as f:
    parameter_values = yaml.safe_load(f)

# Directories
RESULT_DIR = parameter_values["result_dir"]
WORKING_DIR = parameter_values["working_dir"]
FASTQ_DIR = parameter_values["fastq_dir"]
TRIMMED_DIR = WORKING_DIR + "trimmed/"

# Create directories
os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(WORKING_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)

# FASTP parameters
TRIMMING_THRESHOLD = str(parameter_values["qualified_quality_phred"]) # convert to string for subprocess call bash command 

# MUTMAP parameters
REF_GENOME_FASTA = parameter_values["ref_genome"]

# General parameters
THREADS = str(parameter_values["threads"])

#{'result_dir': 'results/', 'working_dir': 'temp/', 'units': 'config/samples.tsv', 'refs': {'genome': 'config/refs/mutmap_ref.fasta'}, 'fastp': {'qualified_quality_phred': 30}, 'mutmap': {'window_size': '2000', 'step_size': '100', 'max_depth': '250', 'min-depth': '8', 'n-rep': '5000', 'min-mq': '40', 'min-bq': '18', 'n_ind': '20'}, 'snpeff': {'database_url': 'https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_Arabidopsis_thaliana.zip', 'database': 'Arabidopsis_thaliana'}}
#print(parameter_values["result_dir"])

########################
# Import list of samples
########################

samples_df = pd.read_csv(parameter_values["samples"], index_col = "sample")
samples = samples_df.index.to_list()

###############
# Read trimming
###############

for sample_name in samples:
    print("################################")
    print("Trimming reads of:", sample_name)
    print("#################################")
    original_fastq_forward = samples_df.loc[sample_name, "fq1"]
    original_fastq_reverse = samples_df.loc[sample_name, "fq2"]
    if sample_is_single_end(sample_name): # single-end >> empty fq2 field 
        bash_trimming_cmd = "fastp -q " + TRIMMING_THRESHOLD + " --thread " + THREADS + " -i " + original_fastq_forward + " -o " + TRIMMED_DIR + sample_name + "_trimmed_R1.fq.gz" 
        subprocess.call(bash_trimming_cmd, shell=True)
    else:
        bash_trimming_cmd = "fastp -q " + TRIMMING_THRESHOLD + " --thread " + THREADS + " -i " + original_fastq_forward + " -I " + original_fastq_reverse + " -o " + TRIMMED_DIR + sample_name + "_trimmed_R1.fq.gz" + " -O " + TRIMMED_DIR + sample_name + "_trimmed_R2.fq.gz"        
        subprocess.call(bash_trimming_cmd, shell=True)

#########
## MutMap
#########
testfiles = get_trimmed_files_per_sample_type(trimmed_dir=TRIMMED_DIR, sample_type="cultivar")


for sample_name in samples:
    # Get the fastq files of the 'cultivar' (the reference genotype)
    cultivar_files = get_trimmed_files_per_sample_type(TRIMMED_DIR, sample_type="cultivar")
    mutant_files = get_trimmed_files_per_sample_type(TRIMMED_DIR, sample_type="bulk")
    print(cultivar_files)
    print(mutant_files)

    print("################################")
    print("MutMap analysis of ", sample_name)
    print("#################################")
    """   if sample_is_single_end(sample_name):
            trimmed_fastq_forward = trimmed_file_dir + sample_name + "_trimmed_R1.fq.gz"
            bash_cmd = "mutmap -r " + REF_GENOME_FASTA + " -c " + trimmed_fastq_forward """

############
# Clean up
#############
cleaning_cmd = "rm -r " + WORKING_DIR
#subprocess.call(cleaning_cmd, shell=True)