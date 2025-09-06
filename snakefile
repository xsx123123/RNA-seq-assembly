# Author : jian zhang 
# *---utf-8---*
# Date   : 2025-8-19
# Version: 1.0
# --------- main snakefile --------- #
configfile: "config/config.yaml"
configfile: "config.yaml"

# include all rules from the rules directory
include: 'rules/log.smk'
include: 'rules/id_convert.smk'
# include: 'rules/qc.smk'
# include: 'rules/trim.smk'
# include: 'rules/long-read-qc.smk'
# include: 'rules/STAR.smk'
# include: 'rules/rnabloom.smk'
# include: 'rules/Corset.smk'
include: 'rules/salmon.smk'

# Set the default target to run
rule all:
    input:
        # Short-Read raw data QC and trimming reports
        # "../01.qc/multiqc/",
        # "../01.qc/multiqc_trim/",
        # expand("../01.qc/trim/{sample}.R1.fastp.fq.gz", sample=load_samples.keys()),
        # expand("../01.qc/trim/{sample}.R2.fastp.fq.gz", sample=load_samples.keys()),
        # expand("../01.qc/trim/{sample}.fastp.html", sample=load_samples.keys()),
        # expand("../01.qc/trim/{sample}.fastp.json", sample=load_samples.keys()),
        # Long-Read raw data QC and trißmming reports
        # expand("../01.qc/nanoplot/{long_sample}/NanoPlot-report.html",long_sample=long_read_samples.keys()),
        # expand("../01.qc/long_read_trim/{long_sample}.fastplong.fq.gz",long_sample=long_read_samples.keys()),
        # expand("../01.qc/long_read_trim/{long_sample}.fastplong.html",long_sample=long_read_samples.keys()),
        # expand("../01.qc/long_read_trim/{long_sample}.fastplong.json",long_sample=long_read_samples.keys()),
        # STAR mapping and BAM index
        # expand("../02.star/{sample}/{sample}.sort.bam",sample=load_samples.keys()),
        # expand("../02.star/{sample}/{sample}.sort.bam.bai",sample=load_samples.keys()),
        # "../02.star/multiqc_star/",
        # RNA-Bloom assembly
        # "../03.assembly/rnabloom_long_assembly"
        # expand("../03.bowtie/{sample}/{sample}.bam",sample=load_samples.keys()),
        # expand("../04.expression/{sample}.corset-reads",sample=load_samples.keys()),
        expand("../04.mapping/salmon_quant/{sample}/",sample=load_samples.keys()),
# --------- main snakefile --------- #