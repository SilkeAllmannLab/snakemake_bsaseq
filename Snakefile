#########################################
# Snakemake pipeline for RNA-Seq analysis
#########################################


###########
# Libraries
###########
import pandas as pd
import os
from subprocess import check_output

###############
# Configuration
###############

configfile: 'config/config.yaml' # where to find parameters

# Directories
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]
TRIMMED_DIR = WORKING_DIR + "trimmed/"

# MutMap parameters
REF_GENOME = config["ref_genome"]

########################
# Samples and conditions
########################

# create lists containing the sample names and conditions
samples = pd.read_csv(config["samples"], dtype=str, index_col=0, sep=",")
SAMPLES = samples.index.values.tolist()
REF_SAMPLE_NAME = samples.query("sample_type == 'reference'").index.values.tolist()
MUTANTS = [s for s in SAMPLES if s not in REF_SAMPLE_NAME]

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

def get_variant_file_of_reference_sample():
    """
    This function re-creates the path to the VCF file of the reference sample.
 
    Returns
    -------
    A string with the path to the reference VCF file 
    """
    # Verify that 'reference' is in the 'sample_type' column of the Pandas samples_df
    sample_types = samples.sample_type.unique()
    #if "reference" not in sample_types:
    #    raise ValueError("sample type should be either equal to 'reference' or 'mutant'")
    
    # Get the sample name corresponding to the reference sample
    ref_sample_name = samples.query("sample_type == 'reference'").index.values[0]

    # Re-create the path of the VCF file
    vcf_file = WORKING_DIR + "vcf/" + ref_sample_name + ".vcf.gz"
    return vcf_file

def get_number_of_individuals_of_sample(wildcards):
    """
    This function get the number of individuals for a given sample.
 
    Returns
    -------
    A string with the number (integer) of individuals per sample 
    """
    samples_filtered = samples.filter(items = [wildcards], axis=0)
    print(samples_filtered)
    nb_of_individuals = samples_filtered.n_individuals
    return int(nb_of_individuals) # conversion to integer for compatibility with MutMap

#################
# Desired outputs
#################
MULTIQC = RESULT_DIR + "multiqc_report.html"
#BAMS = expand(WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.coord_sorted.dedup.bam", sample=SAMPLES)
MUTMAP_VCF = expand(RESULT_DIR + "{sample}/mutmap/30_vcf/mutmap.vcf.gz", sample=SAMPLES)
MUTMAP_ANNOTATED_VCF = expand(RESULT_DIR + "{sample}/snpeff/mutmap_annotated.vcf.gz", sample=SAMPLES)

VCF = expand(RESULT_DIR + "merged_vcf/{sample}.merged_with_reference.vcf", sample=SAMPLES)
MUTPLOTS = expand(RESULT_DIR + "mutplot/{mutant}/mutmap_plot.png", mutant=MUTANTS)

rule all:
    input:
        MULTIQC,
        VCF,
        MUTPLOTS
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
        WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.bam"
    output:
        temp(WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.coord_sorted.bam")
    message:
        "sorting {wildcards.sample} bam file by coordinate"
    threads: 4
    shell:
        "samtools sort -@ {threads} {input} > {output}"

rule mark_duplicate:
    input:
        WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.coord_sorted.bam"
    output:
        WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.coord_sorted.dedup.bam"
    message:
        "marking duplicates in {wildcards.sample} bam file"
    threads: 4
    shell:
        "samtools markdup -@ {threads} {input} {output}"


##############
# Variant call
##############

rule call_variants:
    input:
        bam = WORKING_DIR + "mapped/{sample}.qname_sorted.fixed.coord_sorted.bam"
    output:
        vcf = WORKING_DIR + "vcf/{sample}.vcf.gz"
    message:
        "Call variants for {wildcards.sample}"
    threads: 5
    params:
        mapping_quality = config["mutmap"]["min_mq"],
        minimum_base_quality = config["mutmap"]["min_bq"],
        adjust_mapping_quality = config["mutmap"]["adjust_mq"],
        ref_genome = REF_GENOME
    shell:
        "bcftools mpileup "
        "--annotate AD,ADF,ADR "                       # allelic depth, allelic depth forward strand, allelic depth reverse strand
        "--no-BAQ "                                    # Disable probabilistic realignment for the computation of base alignment quality (BAQ). 
        "--min-MQ {params.mapping_quality} "           # Minimum mapping quality for an alignment to be used [0]
        "--min-BQ {params.minimum_base_quality} "      # Minimum base quality for a base to be considered [13]
        "--adjust-MQ {params.adjust_mapping_quality} " # "adjust-MQ" in mpileup.
        "--output-type u "                             # uncompressed output
        "--fasta-ref {params.ref_genome} "
        "{input.bam} | "
        "bcftools call -vm -f GQ,GP -O u | "
        "bcftools filter -O z -o {output.vcf}"

rule index_variant_files:
    input:
        vcf = WORKING_DIR + "vcf/{sample}.vcf.gz"
    output:
        index = WORKING_DIR + "vcf/{sample}.vcf.gz.csi"
    message:
        "Indexing {wildcards.sample} VCF file"
    shell:
        "bcftools index {input.vcf}"


rule add_reference_to_variant_file_for_mutplot:
    input:
        vcf = WORKING_DIR + "vcf/{sample}.vcf.gz",
        index = WORKING_DIR + "vcf/{sample}.vcf.gz.csi"
    output:
        vcf_merged = RESULT_DIR + "merged_vcf/{sample}.merged_with_reference.vcf"
    message:
        "Adding reference SNPs to {wildcards.sample} vcf file"
    params:
        reference_vcf = get_variant_file_of_reference_sample()
    shell:
        "bcftools merge "
        "--force-samples "          # for the reference sample that will be duplicated
        "--output-type v "          # uncompressed vcf
        "{params.reference_vcf} "   # name of the reference sample
        "{input.vcf} > {output.vcf_merged}"
    
########
# MutMap
########

rule mutplot:
    input:
        ref_fasta = REF_GENOME,
        vcf = RESULT_DIR + "merged_vcf/{mutant}.merged_with_reference.vcf"
    output:
        mutplot = RESULT_DIR + "mutplot/{mutant}/mutmap_plot.png"
    message:
        "Running mutplot on {wildcards.mutant}"
    params:
        window_size           = config["mutmap"]["window_size"],
        step_size             = config["mutmap"]["step_size"],
        n_individuals         = lambda wildcards: int(samples.filter(items = [wildcards.mutant], axis=0).iloc[0]["n_individuals"]),
        outdir_mutmap         = RESULT_DIR + "mutplot/{mutant}/"  
    threads: 10
    shell:
        "mutplot "
        "--vcf {input.vcf} "
        "--window {params.window_size} "
        "--N-bulk {params.n_individuals} "
        "--out {params.outdir_mutmap}"

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


###########
# QTLseqR
###########
# picard CreateSequenceDictionary -R config/refs/mutmap_ref.fasta -O config/refs/mutmap_ref.dict
# samtools faidx config/refs/mutmap_ref.fasta

#gatk3 -T VariantsToTable -V results/merged_vcf/mutant_01.merged_with_reference.vcf 
# -F CHROM -F POS -F REF -F ALT -GF AD -GF DP -GF GQ -GF PL -R 
# config/refs/mutmap_ref.fasta -o results/mutant_01.table


############################
# Plots of mapping statistics
#############################




