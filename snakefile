# Author : jian zhang 
# *---utf-8---*
# Date   : 2025-8-19
# Version: 1.0

configfile: "config/config.yaml"
configfile: "user_config.yaml"

include: 'rules/ID_Convert.smk'
include: 'rules/qc.smk'
include: 'rules/trim.smk'

rule all:
    input:
        "../01.qc/multiqc/multiqc_report.html",
        expand("../01.qc/trim/{sample}.R1.fastp.fq.gz", sample=load_samples.keys()),
        expand("../01.qc/trim/{sample}.R2.fastp.fq.gz", sample=load_samples.keys()),
        expand("../01.qc/trim/{sample}.fastp.html", sample=load_samples.keys()),
        expand("../01.qc/trim/{sample}.fastp.json", sample=load_samples.keys()),