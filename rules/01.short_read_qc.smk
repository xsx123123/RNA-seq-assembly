#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
# logger.info('Running FastQC for quality control on R1 raw sequencing data')
rule short_read_qc_r1:
    input:
        r1 = os.path.join(config["raw_data_path"], "{sample}.R1.fq.gz"),
    output:
        r1_html = "../01.qc/short_read_qc_r1/{sample}.R1_fastqc.html",
        r1_zip = "../01.qc/short_read_qc_r1/{sample}.R1_fastqc.zip", 
    conda:
        "../envs/fastqc.yaml",
    log:
        r1 = "../logs/short_read_qc_r1/{sample}.r1.fastqc.log",
    benchmark:
        "../benchmarks/{sample}.r1.fastqc.txt",
    params:
        out_dir = "../01.qc/short_read_qc_r1/",
    message:
        "Running FastQC on {input.r1}",
    threads: 1
    shell:
        """
        fastqc {input.r1} \
               -o {params.out_dir} \
               --threads {threads} &> {log.r1}
        """

# logger.info('Running FastQC for quality control on R2 raw sequencing data')
rule short_read_qc_r2:
    input:
        r2 = os.path.join(config["raw_data_path"], "{sample}.R2.fq.gz"),
    output:
        r2_html = "../01.qc/short_read_qc_r2/{sample}.R2_fastqc.html",
        r2_zip = "../01.qc/short_read_qc_r2/{sample}.R2_fastqc.zip", 
    conda:
        "../envs/fastqc.yaml",
    log:
        r2 = "../logs/short_read_qc_r2/{sample}.r2.fastqc.log",
    benchmark:
        "../benchmarks/{sample}.r2.fastqc.txt",
    params:
        out_dir = "../01.qc/short_read_qc_r2",
    message:
        "Running FastQC on {input.r2}",
    threads: 1
    shell:
        """
        fastqc {input.r2} \
               -o {params.out_dir} \
               --threads {threads} &> {log.r2}
        """

# logger.info('Run MultiQC to summarize R1 fastqc QC reports')
rule short_read_multiqc_r1:
    input:
        fastqc_files_r1 = expand("../01.qc/short_read_qc_r1/{sample}.R1_fastqc.zip", sample=load_samples.keys()),
    output:
        report_dir = directory("../01.qc/short_read_r1_multiqc/")
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to aggregate R1 FastQC reports",
    params:
        fastqc_reports = "../01.qc/short_read_qc_r1",
        report = "multiqc_r1_raw-data_report.html",
        title = "r1-raw-data-multiqc-report",
    log:
        "../logs/qc/multiqc-r1.log",
    benchmark:
        "../benchmarks/multiqc-r1.txt",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """

# logger.info('Run MultiQC to summarize R2 fastqc QC reports')
rule short_read_multiqc_r2:
    input:
        fastqc_files_r2 = expand("../01.qc/short_read_qc_r2/{sample}.R2_fastqc.zip", sample=load_samples.keys()),
    output:
        report_dir = directory("../01.qc/short_read_r2_multiqc/")
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to aggregate R2 FastQC reports",
    params:
        fastqc_reports = "../01.qc/short_read_qc_r2",
        report = "multiqc_r2_raw-data_report.html",
        title = "r2-raw-data-multiqc-report",
    log:
        "../logs/qc/multiqc-r2.log",
    benchmark:
        "../benchmarks/multiqc-r2.txt",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """
# ----- rule ----- #