#!/usr/bin/env python3

from random import sample
import yaml
import pandas as pd
import os
import subprocess

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
MUTMAP_DIR = RESULT_DIR + "mutmap/"

# Clean-up existing directories and create fresh ones
subprocess.call("rm -r " + RESULT_DIR, shell=True)
subprocess.call("rm -r " + WORKING_DIR, shell=True)
os.makedirs(RESULT_DIR, exist_ok=True)
os.makedirs(WORKING_DIR, exist_ok=True)
os.makedirs(TRIMMED_DIR, exist_ok=True)
os.makedirs(RESULT_DIR + "fastp/", exist_ok=True)

# FASTP parameters
TRIMMING_THRESHOLD = str(parameter_values["qualified_quality_phred"]) # convert to string for subprocess call bash command 

# MUTMAP parameters
REF_GENOME_FASTA = parameter_values["ref_genome"]
MEM_MUTMAP = "1G"

# SNPEFF parameters
SNPEFF_FORMAT = parameter_values["snpeff"]["output_format"]
SNPEFF_DB = parameter_values["snpeff"]["database"] 

# General parameters
THREADS = str(parameter_values["threads"])

#{'result_dir': 'results/', 'working_dir': 'temp/', 'units': 'config/samples.tsv', 'refs': {'genome': 'config/refs/mutmap_ref.fasta'}, 'fastp': {'qualified_quality_phred': 30}, 'mutmap': {'window_size': '2000', 'step_size': '100', 'max_depth': '250', 'min-depth': '8', 'n-rep': '5000', 'min-mq': '40', 'min-bq': '18', 'n_ind': '20'}, 'snpeff': {'database_url': 'https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_Arabidopsis_thaliana.zip', 'database': 'Arabidopsis_thaliana'}}
#print(parameter_values["result_dir"])

########################
# Import list of samples
########################

samples_df = pd.read_csv(parameter_values["samples"], index_col = "sample_name")
samples = samples_df.index.to_list()


##################
# Helper functions
##################
def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples_df.columns:
        return True
    else:
        return pd.isnull(samples_df.loc[(sample), "fq2"])

def get_trimmed_files_of_reference_sample():
    """
    This function re-creates the path to the trimmed files of the reference sample.
 
    Returns
    -------
    A string with the ',' and the reference trimmed file path (input for MutMap)
    """
    # Verify that 'reference' is in the 'sample_type' column of the Pandas samples_df
    sample_types = samples_df.sample_type.unique()
    if "reference" not in sample_types:
        raise ValueError("sample type should be either equal to 'reference' or 'mutant'")
    
    # Get the sample name corresponding to the reference sample
    ref_sample_name = samples_df.query("sample_type == 'reference'").index.values[0]

    # Re-create the path of the trimmed file
    if sample_is_single_end(ref_sample_name):
        trimmed_files = ",".join(TRIMMED_DIR + ref_sample_name + "_trimmed_R1.fq")
    else:
        trimmed_files = ",".join([TRIMMED_DIR + ref_sample_name + "_trimmed_R1.fq", TRIMMED_DIR + ref_sample_name + "_trimmed_R2.fq"]) 
    
    return trimmed_files

def get_trimmed_files_of_sample(sample):
    """
    This function:
      1. Re-creates the files names of the trimmed files for a given mutant file 
      2. Returns a list with one or two elements. 

    Parameters
    ---------
    sample: str
      name of the sample as given in the sample_name column of the 'samples.csv' file

    Returns
    -------
    A string of the path to trimmed file(s) names with a comma ',' as separator for the files (input for MutMap)
    """
    # Verify that the sample have 'mutant' specified as sample_type
    if samples_df.loc[sample_name, "sample_type"] != "mutant":
        raise ValueError("Please verify that sample ", sample, " has 'mutant' specified as 'sample_type")
    else:
        # Re-create the path of the trimmed file
        if sample_is_single_end(sample):
            trimmed_files = ",".join(TRIMMED_DIR + sample + "_trimmed_R1.fq")
        else:
            trimmed_files = ",".join([TRIMMED_DIR + sample + "_trimmed_R1.fq", TRIMMED_DIR + sample + "_trimmed_R2.fq"]) 
    return trimmed_files

def get_number_of_individuals(sample):
    """
    Imports the samples .csv dataframe Returns the number of mutant individuals used for sequencing
    
    Parameter
    ---------
    sample: str
        The name of the mutant used to return its corresponding number of individuals
    
    Returns
    -------
    An integer corresponding to the number of individuals
    """
    if sample not in samples:
        raise NameError(sample," is not an element of your list of samples (sample_name column in .csv)")

    n_ind = int(samples_df.loc[sample, "n_individuals"])
    return n_ind 

################
# Pipeline steps
################
        
###################
# 01. Read trimming
###################

for sample_name in samples:
    print("################################")
    print("Trimming reads of:", sample_name)
    print("#################################")
    original_fastq_forward = samples_df.loc[sample_name, "fq1"]
    original_fastq_reverse = samples_df.loc[sample_name, "fq2"]
    if sample_is_single_end(sample_name): # single-end >> empty fq2 field 
        bash_trimming_cmd = "fastp -q " + TRIMMING_THRESHOLD + " --thread " + THREADS + " -i " + original_fastq_forward + " -o " + TRIMMED_DIR + sample_name + "_trimmed_R1.fq" 
        subprocess.call(bash_trimming_cmd, shell=True)
    else:
        bash_trimming_cmd = "fastp -q " + TRIMMING_THRESHOLD + " --thread " + THREADS + " -i " + original_fastq_forward + " -I " + original_fastq_reverse + " -o " + TRIMMED_DIR + sample_name + "_trimmed_R1.fq" + " -O " + TRIMMED_DIR + sample_name + "_trimmed_R2.fq"        
        subprocess.call(bash_trimming_cmd, shell=True)

    subprocess.call("mv fastp.html " + RESULT_DIR + "fastp/" + sample_name + ".html", shell=True)
    subprocess.call("mv fastp.json " + RESULT_DIR + "fastp/" + sample_name + ".json", shell=True)

#############
## 02. MutMap
#############

# Sample name of reference sample
ref_sample_name = samples_df.query("sample_type == 'reference'").index.values[0]

for sample in samples:
    if sample != ref_sample_name:
        print("########################################")
        print("Performing mutmap analysis for:", sample)
        print("########################################")
        ref_trimmed_fastq_files = get_trimmed_files_of_reference_sample()
        mutant_trimmed_fastq_files = get_trimmed_files_of_sample(sample=sample)
        n_ind = get_number_of_individuals(sample=sample)

        mutmap_cmd = "mutmap --ref " + REF_GENOME_FASTA +  " -c " + ref_trimmed_fastq_files + " -b " + mutant_trimmed_fastq_files + " -n " + str(n_ind) + " -o " + "mutmap/" + " -t " + THREADS + " --mem " + MEM_MUTMAP
        subprocess.call(mutmap_cmd, shell=True)
        subprocess.call("mv mutmap/" + " " + RESULT_DIR + sample, shell=True) 


############
# 03. snpEff
############

for sample in samples:
    if sample != ref_sample_name:
        print("##############################")
        print("Annotating SNPs for:", sample)
        print("##############################")
        
        # decompress variant file
        vcf_decompress_cmd = "gzip -d " + RESULT_DIR + sample + "/30_vcf/mutmap.vcf.gz"
        subprocess.call(vcf_decompress_cmd, shell=True)

        # Annotate SNPs using snpEff
        snpeff_cmd = "snpEff Arabidopsis_thaliana -stats " + RESULT_DIR + sample + "/snpeff_stats.html " + RESULT_DIR + sample + "/30_vcf/mutmap.vcf > " + RESULT_DIR + sample + "/mutmap_annotated.vcf"   
        print(snpeff_cmd)
        subprocess.call(snpeff_cmd, shell=True)    

#################
# End of pipeline
#################