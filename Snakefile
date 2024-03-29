#########################################
# Snakemake pipeline for BSA-seq analysis
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

# Reference genome
REF_GENOME = config["ref_genome"]

########################
# Samples and conditions
########################

# create lists containing the sample names and conditions
samples = pd.read_csv(config["samples"], dtype=str, index_col=0, sep=",")
SAMPLES = samples.index.values.tolist()


####################################################################
# Wildcards constraint (to deal with issues with dots in file names)
####################################################################

#wildcard_constraints:
#    sample="[0-9][Aa-Zz]+" # dots are not allowed > https://stackoverflow.com/questions/4742357/how-to-match-only-strings-that-do-not-contain-a-dot-using-regular-expressions

###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    "This function detect missing value in the column 2 of the units.tsv"
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

#################
# Desired outputs
#################
MULTIQC = RESULT_DIR + "multiqc_report.html"
BAMS = expand(WORKING_DIR + "mapped/{sample}.bam", sample=SAMPLES)
CLEAN_BAMS = expand(WORKING_DIR + "samtools/dedup/{sample}_qname_sorted_fixed_coord_sorted_dedup.bam", sample=SAMPLES)

GCVFs = expand(WORKING_DIR + "gatk/{sample}.g.vcf.gz", sample=SAMPLES)
VCF = RESULT_DIR + "vcf/all_samples.vcf.gz"
VCF_TABLE = RESULT_DIR + "all_samples.variants.tsv"

rule all:
    input:
         MULTIQC, BAMS, CLEAN_BAMS,
         GCVFs, 
         VCF, 
         VCF_TABLE # suitable for QTLseqR
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
        WORKING_DIR + "fastp/{sample}.log.txt"
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
        if sample_is_single_end(wildcards.sample):
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
        WORKING_DIR + "mapped/{sample}.bam"
    message:
        "mapping {wildcards.sample} reads to genomic reference"
    params:
        db_prefix = WORKING_DIR + "index/genome"
    threads: 10
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
            shell("bwa mem -v 0 -t {threads} -R '@RG\\tID:{READ_GROUP}\\tPL:ILLUMINA\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}' {params.db_prefix} {input.forward_fastq} >{output}")
        else:
            shell("bwa mem -v 0 -t {threads} -R '@RG\\tID:{READ_GROUP}\\tPL:ILLUMINA\\tLB:{wildcards.sample}\\tSM:{wildcards.sample}' {params.db_prefix} {input.forward_fastq} {input.reverse_fastq} >{output}")

###############################################
# Post-alignment steps prior to variant calling
###############################################

rule samtools_sort_by_qname:
    input:
        WORKING_DIR + "mapped/{sample}.bam"
    output:
        WORKING_DIR + "samtools/sort_qname/{sample}_qname_sorted.bam"
    message:
         "sorting {wildcards.sample} bam file by read name (QNAME field)"
    threads: 10
    shell:
        "samtools sort -n -@ {threads} {input} > {output}"

rule samtools_fixmate:
    input:
        WORKING_DIR + "samtools/sort_qname/{sample}_qname_sorted.bam"
    output:
        WORKING_DIR + "samtools/fixmate/{sample}_qname_sorted_fixed.bam"
    message:
        "Fixing mate in {wildcards.sample} sorted bam file"
    threads: 10
    shell:
        "samtools fixmate -m -@ {threads} {input} {output}"

rule samtools_sort_by_coordinates:
    input:
        WORKING_DIR + "samtools/fixmate/{sample}_qname_sorted_fixed.bam"
    output:
        WORKING_DIR + "samtools/sort_coords/{sample}_qname_sorted_fixed_coord_sorted.bam"
    message:
        "sorting {wildcards.sample} bam file by coordinate"
    threads: 10
    shell:
        "samtools sort -@ {threads} {input} > {output}"


rule mark_duplicate:
    input:
        WORKING_DIR + "samtools/sort_coords/{sample}_qname_sorted_fixed_coord_sorted.bam"
    output:
        WORKING_DIR + "samtools/dedup/{sample}_qname_sorted_fixed_coord_sorted_dedup.bam"
    message:
        "marking duplicates in {wildcards.sample} bam file"
    threads: 10
    shell:
        "samtools markdup -@ {threads} {input} {output};"
        "samtools index {output}"

########################
# Variant call with GATK
########################
rule prepare_fasta_for_gatk:
    input:
        ref = REF_GENOME
    output:
        ref_dict = os.path.splitext(REF_GENOME)[0] + ".dict"
    message:
        "Creating sequence dictionary and index for {REF_GENOME}"
    shell:
        "samtools faidx {input.ref};"
        "picard CreateSequenceDictionary -R {input.ref} -O {output.ref_dict}"

rule call_variants_with_gatk:
    input:
        bam = WORKING_DIR + "samtools/dedup/{sample}_qname_sorted_fixed_coord_sorted_dedup.bam",
        ref = REF_GENOME,
        ref_dict = rules.prepare_fasta_for_gatk.output.ref_dict
    output:
        gvcf = WORKING_DIR + "gatk/{sample}.g.vcf.gz"
    message:
        "Calling {wildcards.sample} variants (GATK HaplotypeCaller)"
    shell:
        "gatk HaplotypeCaller  "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.gvcf} "
        "-ERC GVCF "
 
rule create_list_of_gvcf_files_for_joint_genotyping:
    input:
        gvcfs = expand(WORKING_DIR + "gatk/{sample}.g.vcf.gz", sample=SAMPLES) 
    output:
        gvcf_list = WORKING_DIR + "gatk/gvcfs.list"
    params:
        path_to_gvcf_files = WORKING_DIR + "gatk/"
    message:
        "Creating list of individual GVCF files for joint genotyping"
    shell:
        "find {params.path_to_gvcf_files} -type f -name '*.g.vcf.gz' > {output}"


rule combine_gvcf_files:
    input:
        gvcf_list = WORKING_DIR + "gatk/gvcfs.list",
        ref =       REF_GENOME
    output:
        combined_gvcf = WORKING_DIR + "gatk/all_samples.gatk.vcf.gz"
    message:
        "Combine individual GVCF files (GATK CombineGVCFs)"
    shell:
        "gatk CombineGVCFs "
        "-R {input.ref} "
        "--variant {input.gvcf_list} "
        "-O {output.combined_gvcf}"

rule joint_genotypying_with_gatk:
    input:
        combined_gvcf = WORKING_DIR + "gatk/all_samples.gatk.vcf.gz",
        ref =  REF_GENOME
    output:
        vcf = RESULT_DIR + "vcf/all_samples.vcf.gz"
    message:
        "Joint cohort variant calling (GATK GenotypeGVCFs)"
    shell:
        "gatk GenotypeGVCFs "
        "-R {input.ref} "
        "--variant {input.combined_gvcf} "
        "-O {output.vcf}"


###########################
# prepare table for QTLseqR
###########################

rule convert_variants_to_table:
    input:
        picard_dict = rules.prepare_fasta_for_gatk.output.ref_dict, # picard dict
        vcf = RESULT_DIR + "vcf/all_samples.vcf.gz",
        ref_genome = REF_GENOME
    output:
        table = RESULT_DIR + "all_samples.variants.tsv"
    message:
        "Converting VCF from all samples to table for QTLseqR"
    shell:
        "gatk VariantsToTable "
        "-V {input.vcf} "
        "-F CHROM -F POS -F REF -F ALT "
        "-GF AD -GF DP -GF GQ -GF PL "
        "-R {input.ref_genome} "
        "-O {output.table}"