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


# Clean potential left-over mutmap folders
if os.path.isdir(WORKING_DIR + "mutmap"):
    os.rmdir(WORKING_DIR + "mutmap")


########################
# Samples and conditions
########################

# create lists containing the sample names and conditions
samples = pd.read_csv(config["samples"], dtype=str, index_col=0, sep=",")
SAMPLES = samples.index.values.tolist()


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

def get_trimmed_files_per_sample_type(trimmed_dir = WORKING_DIR + "trimmed", sample_type = "cultivar"):
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

#################
# Desired outputs
#################
MULTIQC = RESULT_DIR + "multiqc_report.html"
MUTMAP_VCF = RESULT_DIR + "mutmap/30_vcf/mutmap.vcf.gz"
MUTMAP_ANNOTATED_VCF = RESULT_DIR + "snpeff/mutmap_annotated.vcf.gz"

rule all:
    input:
        MULTIQC,
        MUTMAP_VCF,
        MUTMAP_ANNOTATED_VCF
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
        ref_fasta = config["ref_genome"],
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
        "mutmap --ref {input.ref_fasta} "
        "--cultivar {params.cultivar_files} "
        "--bulk {params.bulk_files} "          
        "--threads {threads}  "
        "--window {params.window_size} "
        "--N-bulk {params.n_individuals} "
        "--out {params.outdir_mutmap};"
        "cp -r {params.outdir_mutmap} {params.outdir_final}"


########
# snpEff
########

rule snpeff:
    input:
        RESULT_DIR + "mutmap/30_vcf/mutmap.vcf.gz"
    output:
        csv = RESULT_DIR + "snpeff/stats.csv",
        ann = RESULT_DIR + "snpeff/mutmap_annotated.vcf.gz"
    message:
        "Annotating SNPs with snpEff"
    params:
        snpeff_db = config["snpeff"]["database"],
        output_format = config["snpeff"]["output_format"]
    shell:
        "snpEff ann "
        "-o {params.output_format} "
        "-csvStats "                            # creates CSV summary file instead of HTML
        "-stats {output.csv} "
        "{params.snpeff_db} "
        "{input} > {output.ann}"

############################
# Plots of mapping statistics
#############################




