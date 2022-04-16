#########################################
# Snakemake pipeline for RNA-Seq analysis
#########################################


###########
# Libraries
###########
from random import sample
import pandas as pd
import os
from subprocess import check_output

###############
# Configuration
###############

configfile: "config/config.yaml" # where to find parameters

# Directories
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]
TRIMMED_DIR = WORKING_DIR + "trimmed/"

# MutMap parameters
REF_GENOME = config["ref_genome"]

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
    if "fq2" in samples.columns:
        return False
    else:
        return False
        #return pd.isnull(samples.loc[(sample), "fq2"])

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

def get_trimmed_files_of_sample(wildcards):
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
    if sample_is_single_end(wildcards):
        trimmed_files = ",".join(TRIMMED_DIR + wildcards.sample + "_trimmed_R1.fq")
    else:
        trimmed_files = ",".join([TRIMMED_DIR + wildcards.sample + "_trimmed_R1.fq", TRIMMED_DIR + wildcards.sample + "_trimmed_R2.fq"]) 
    return trimmed_files

def get_trimmed_files_of_reference_sample():
    """
    This function re-creates the path to the trimmed files of the reference sample.
 
    Returns
    -------
    A string with the ',' and the reference trimmed file path (input for MutMap)
    """
    # Verify that 'reference' is in the 'sample_type' column of the Pandas samples_df
    sample_types = samples.sample_type.unique()
    #if "reference" not in sample_types:
    #    raise ValueError("sample type should be either equal to 'reference' or 'mutant'")
    
    # Get the sample name corresponding to the reference sample
    ref_sample_name = samples.query("sample_type == 'reference'").index.values[0]

    # Re-create the path of the trimmed file
    if sample_is_single_end(ref_sample_name):
        trimmed_files = ",".join(TRIMMED_DIR + ref_sample_name + "_trimmed_R1.fq")
    else:
        trimmed_files = ",".join([TRIMMED_DIR + ref_sample_name + "_trimmed_R1.fq", TRIMMED_DIR + ref_sample_name + "_trimmed_R2.fq"]) 
    
    return trimmed_files



#################
# Desired outputs
#################
MULTIQC = RESULT_DIR + "multiqc_report.html"
BAMS = expand(WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.coord_sorted.dedup.bam", sample=SAMPLES)
MUTMAP_VCF = expand(RESULT_DIR + "{sample}/mutmap/30_vcf/mutmap.vcf.gz", sample=SAMPLES)
MUTMAP_ANNOTATED_VCF = expand(RESULT_DIR + "{sample}/snpeff/mutmap_annotated.vcf.gz", sample=SAMPLES)


rule all:
    input:
        MULTIQC,
        BAMS
#        MUTMAP_VCF,
#        MUTMAP_ANNOTATED_VCF
    message:
        "BSA-seq pipeline run complete!"
    shell:
        "cp config/config.yaml {RESULT_DIR};"
        "cp config/samples.csv {RESULT_DIR};"
        #"rm -r {WORKING_DIR}"

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


###########
# Alignment
###########

rule bwa_index:
    input:
        genome = REF_GENOME
    output:
        WORKING_DIR + "index/genome.amb",
        WORKING_DIR + "index/genome.ann",
        WORKING_DIR + "index/genome.pac",
        WORKING_DIR + "index/genome.sa",
        WORKING_DIR + "index/genome.bwt" 
    message:
        "building BWA index for the genomic reference"
    params:
        db_prefix = WORKING_DIR + "index/genome"
    shell:
        "bwa index -p {params.db_prefix} {input}"

rule uncompress:
    input:
        forward_fastq = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq.gz",
        reverse_fastq = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq.gz"
    output:
        forward_fastq = WORKING_DIR + "trimmed/{sample}_forward.fastq",
        reverse_fastq = WORKING_DIR + "trimmed/{sample}_reverse.fastq"
    message:
        "uncompressing {wildcards.sample} reads"
    run:
        if sample_is_single_end("{wildcards.sample}"):
            shell("gzip -cd {input.forward_fastq} > {output.forward_fastq};touch {output.reverse_fastq}")
        else:
            shell("gzip -cd {input.forward_fastq} > {output.forward_fastq}")
            shell("gzip -cd {input.reverse_fastq} > {output.reverse_fastq}")


rule bwa_align:
    input:
        index = [WORKING_DIR + "index/genome." + ext for ext in ["amb","ann","pac","sa","bwt"]],
        forward_fastq = WORKING_DIR + "trimmed/{sample}_forward.fastq",
        reverse_fastq = WORKING_DIR + "trimmed/{sample}_reverse.fastq"
    output:
        temp(WORKING_DIR + "mapped/{sample}.bam")
    message:
        "mapping {wildcards.sample} reads to genomic reference"
    params:
        db_prefix = WORKING_DIR + "index/genome"
    threads: 5
    run:
        # Building the read group id (sequencer_id + flowcell_name + lane_number + barcode)
        SEQUENCER_ID=check_output("head -n 1 " + input.forward_fastq + " |cut -d: -f1",shell=True).decode().strip()
        FLOWCELL_NAME=check_output("head -n 1 " + input.forward_fastq + " |cut -d: -f3",shell=True).decode().strip()
        FLOWCELL_LANE=check_output("head -n 1 " + input.forward_fastq + " |cut -d: -f4",shell=True).decode().strip()
        BARCODE=check_output("head -n 1 " + input.forward_fastq + " |cut -d' ' -f2 |cut -d: -f4",shell=True).decode().strip()
        # Feeding the READ_GROUP_ID to bwa
        READ_GROUP = SEQUENCER_ID + "." + FLOWCELL_NAME + "." + FLOWCELL_LANE + "." + BARCODE
        # If sample is single end, feeding only one fastq file (other outputs an empty BAM file)
        if sample_is_single_end(wildcards.sample):
            shell("bwa mem -v 1 -t {threads} -R '@RG\\tID:{READ_GROUP}\\tPL:ILLUMINA\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}' {params.db_prefix} {input.forward_fastq} >{output}")
        else:
            shell("bwa mem -v 1 -t {threads} -R '@RG\\tID:{READ_GROUP}\\tPL:ILLUMINA\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}' {params.db_prefix} {input.forward_fastq} {input.reverse_fastq} >{output}")

###############################################
# Post-alignment steps prior to variant calling
###############################################

rule samtools_sort_by_qname:
    input:
        WORKING_DIR + "mapped/{sample}.bam"
    output:
        temp(WORKING_DIR + "mapped/{sample}.qname_sorted.bam")
    message:
         "sorting {wildcards.sample} bam file by read name (QNAME field)"
    threads: 4
    shell:
        "samtools sort -n -@ {threads} {input} > {output}"

rule samtools_fixmate:
    input:
        WORKING_DIR + "mapped/{sample}.qname_sorted.bam"
    output:
        temp(WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.bam")
    message:
        "Fixing mate in {wildcards.sample} sorted bam file"
    threads: 4
    shell:
        "samtools fixmate -m -@ {threads} {input} {output}"

rule samtools_sort_by_coordinates:
    input:
        temp(WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.bam")
    output:
        temp(WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.coord_sorted.bam")
    message:
        "sorting {wildcards.sample} bam file by coordinate"
    threads: 4
    shell:
        "samtools sort -@ {threads} {input} > {output}"

rule mark_duplicate:
    input:
        temp(WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.coord_sorted.bam")
    output:
        WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.coord_sorted.dedup.bam"
    message:
        "marking duplicates in {wildcards.sample} bam file"
    threads: 4
    shell:
        "samtools markdup -@ {threads} {input} {output}"

########
# MutMap
########

rule mutmap:
    input:
        ref_fasta = REF_GENOME,
        #mutant_trimmed = get_trimmed_files_of_sample
        #forward_trimmed_files = expand(WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz", sample = SAMPLES),
        #reverse_trimmed_files = expand(WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz", sample = SAMPLES)
    output:
        RESULT_DIR + "{sample}/mutmap/30_vcf/mutmap.vcf.gz"
    message:
        "Running MutMap"
    params:
        cultivar_files        = get_trimmed_files_of_reference_sample(),
        mutant_trimmed = get_trimmed_files_of_sample,
        #bulk_files            = get_trimmed_files_of_sample,
        window_size           = config["mutmap"]["window_size"],
        step_size             = config["mutmap"]["step_size"],
        n_individuals         = config["mutmap"]["n_ind"],
        outdir_mutmap         = "mutmap/",   # only temporary for mutmap
        outdir_final          = RESULT_DIR + "{sample}/mutmap/"
    threads: 10
    shell:
        "mutmap --ref {input.ref_fasta} "
        "--cultivar {params.cultivar_files} "
        "--bulk {params.mutant_trimmed} "          
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
        RESULT_DIR + "{sample}/mutmap/30_vcf/mutmap.vcf.gz"
        #RESULT_DIR + "mutmap/30_vcf/mutmap.vcf.gz"
    output:
        csv = RESULT_DIR + "{sample}/snpeff/stats.csv",
        ann = RESULT_DIR + "{sample}/snpeff/mutmap_annotated.vcf.gz"
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




