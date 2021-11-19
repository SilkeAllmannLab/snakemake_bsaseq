#########################################
# Snakemake pipeline for RNA-Seq analysis
#########################################


###########
# Libraries
###########
import pandas as pd
import os

###############
# Configuration
###############
configfile: "config/config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]


########################
# Samples and conditions
########################

# read the tabulated separated table containing the sample, condition and fastq file information∂DE
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

# create lists containing the sample names and conditions
SAMPLES = units.index.get_level_values('sample').unique().tolist()
samples = pd.read_csv(config["units"], dtype=str,index_col=0,sep="\t")
samplefile = config["units"]


###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
    """This function checks if the sample has paired end or single end reads and returns 1 or 2 names of the fastq files"""
    if sample_is_single_end(wildcards.sample):
        return samples.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trim_names(wildcards):
    """
    This function:
      1. Checks if the sample is paired end or single end
      2. Returns the correct input and output trimmed file names. 
    """
    if sample_is_single_end(wildcards.sample):
        inFile = samples.loc[(wildcards.sample), ["fq1"]].dropna()
        return "--in1 " + inFile[0] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz" 
    else:
        inFile = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
        return "--in1 " + inFile[0] + " --in2 " + inFile[1] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz --out2 "  + WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"

def get_trimmed_files_per_sample_type(trimmed_dir = WORKING_DIR + "trimmed/", sample_type = "cultivar"):
    """
    This function:
      1. Collects the files names of the trimmed files
      2. Filter to keep only the "cultivar" or "bulk" trimmed files.
      3. Returns a list with one or two elements. 
    """
    if not os.path.isdir(trimmed_dir):
        os.mkdir(trimmed_dir)

    trimmed_files = os.listdir(trimmed_dir) # collects all trimmed files
    trimmed_files.sort()                                 # to have R1 before R2

    cultivar_files = [trimmed_dir + f for f in trimmed_files if sample_type in f]
    bulk_files = [trimmed_dir + f for f in trimmed_files if sample_type in f]
    if sample_type == "cultivar":
        return ",".join(cultivar_files) # fastq of cultivar. If you specify fastq, please separate pairs by comma, e.g. -b fastq1,fastq2.
    elif sample_type == "bulk":
        return ",".join(bulk_files)     # fastq of mutant bulk. If you specify fastq, please separate pairs by comma, e.g. -b fastq1,fastq2.
    else:
        "Please check that sample_type == 'cultivar' or sample_type == 'bulk'"

#################
# Desired outputs
#################
MULTIQC = RESULT_DIR + "multiqc_report.html"
MUTMAP_VCF = RESULT_DIR + "mutmap/30_vcf/mutmap.vcf.gz"

rule all:
    input:
        MULTIQC,
        MUTMAP_VCF
    message:
        "BSA-seq pipeline run complete!"
    shell:
        "cp config/config.yaml {RESULT_DIR};"
        "cp config/samples.tsv {RESULT_DIR};"
        "rm -r {WORKING_DIR}"

#######
# Rules
#######


###############
# Read trimming
###############

rule fastp:
    input:
        get_fastq
    output:
        fq1  = temp(WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz"),
        fq2  = temp(WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"),
        html = WORKING_DIR + "fastp/{sample}_fastp.html",
        json = WORKING_DIR + "fastp/{sample}_fastp.json"
    message:"trimming {wildcards.sample} reads"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        in_and_out_files =  get_trim_names,
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    resources: cpus=10
    shell:
        "touch {output.fq2};\
        fastp --thread {threads}  --html {output.html} --json {output.json} \
        --qualified_quality_phred {params.qualified_quality_phred} \
        {params.in_and_out_files} \
        2>{log}"

rule multiqc:
    input:
        expand(WORKING_DIR + "fastp/{sample}_fastp.json", sample = SAMPLES)
    output:
        RESULT_DIR + "multiqc_report.html"
    params:
        fastp_directory = WORKING_DIR + "fastp/",
        outdir = RESULT_DIR
    message: "Summarising fastp reports with multiqc"
    shell:
        "multiqc --force "
        "--outdir {params.outdir} "
        "{params.fastp_directory} "

########
# MutMap
########

rule mutmap:
    input:
        ref_fasta = config["refs"]["genome"],
        forward_trimmed_files = expand(WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz", sample = SAMPLES),
        reverse_trimmed_files = expand(WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz", sample = SAMPLES)
    output:
        RESULT_DIR + "mutmap/30_vcf/mutmap.vcf.gz"
    message:
        "Running MutMap"
    params:
        cultivar_files        = get_trimmed_files_per_sample_type(trimmed_dir = WORKING_DIR + "trimmed/", sample_type = "cultivar"),
        bulk_files            = get_trimmed_files_per_sample_type(trimmed_dir = WORKING_DIR + "trimmed/", sample_type = "bulk"),
        window_size           = config["mutmap"]["window_size"],
        step_size             = config["mutmap"]["step_size"],
        n_individuals         = config["mutmap"]["n_ind"],
        outdir_mutmap         = "mutmap/",   # only temporary for mutmap
        outdir_final          = RESULT_DIR + "mutmap/"
    threads: 10
    shell:
        "rm -r mutmap;"   # remove potential previous mutmap temp directories
        "mutmap --ref {input.ref_fasta} "
        "--cultivar {params.cultivar_files} "
        "--bulk {params.bulk_files} "          
        "--threads {threads}  "
        "--window {params.window_size} "
        "--N-bulk {params.n_individuals} "
        "--out {params.outdir_mutmap};"
        "mv -f {params.outdir_mutmap} {params.outdir_final}"

       

############################
# Plots of mapping statistics
#############################




