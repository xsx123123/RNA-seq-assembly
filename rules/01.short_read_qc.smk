#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
rule short_read_qc_r1:
    input:
        md5_check = "../01.qc/md5_check.tsv",
    output:
        r1_html = "../01.qc/short_read_qc_r1/{sample}_R1_fastqc.html",
        r1_zip = "../01.qc/short_read_qc_r1/{sample}_R1_fastqc.zip",
    conda:
        "../envs/fastqc.yaml",
    log:
        r1 = "../logs/01.short_read_qc_r1/{sample}.r1.fastqc.log",
    params:
        out_dir = "../01.qc/short_read_qc_r1/",
        r1 = os.path.join(config["raw_data_path"],
                          config['convert_md5'],
                          "{sample}",
                          "{sample}" + config['r1_suffix']),
    message:
        "Running FastQC on {wildcards.sample} r1",
    benchmark:
        "../benchmarks/{sample}_r1_fastqc_benchmark.txt",
    threads: 1
    shell:
        """
        fastqc {params.r1} \
               -o {params.out_dir} \
               --threads {threads} &> {log.r1}
        """

rule short_read_qc_r2:
    input:
        md5_check = "../01.qc/md5_check.tsv",
    output:
        r2_html = "../01.qc/short_read_qc_r2/{sample}_R2_fastqc.html",
        r2_zip = "../01.qc/short_read_qc_r2/{sample}_R2_fastqc.zip",
    conda:
        "../envs/fastqc.yaml",
    log:
        r2 = "../logs/01.short_read_qc_r2/{sample}.r2.fastqc.log",
    params:
        out_dir = "../01.qc/short_read_qc_r2",
        r2 = os.path.join(config["raw_data_path"],
                          config['convert_md5'],
                          "{sample}",
                          "{sample}" + config['r2_suffix']),
    message:
        "Running FastQC on {wildcards.sample} r2",
    benchmark:
        "../benchmarks/{sample}_r2_fastqc_benchmark.txt",
    threads: 1
    shell:
        """
        fastqc {params.r2} \
               -o {params.out_dir} \
               --threads {threads} &> {log.r2}
        """

# logger.info('Run MultiQC to summarize R1 fastqc QC reports')
rule short_read_multiqc_r1:
    input:
        fastqc_files_r1 = expand("../01.qc/short_read_qc_r1/{sample}_R1_fastqc.zip", sample=load_samples.keys()),
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
        "../logs/01.multiqc/multiqc-r1.log",
    benchmark:
        "../benchmarks/fastqc_multiqc-r1_benchmark.txt",
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
        fastqc_files_r2 = expand("../01.qc/short_read_qc_r2/{sample}_R2_fastqc.zip", sample=load_samples.keys()),
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
        "../logs/01.multiqc/multiqc-r2.log",
    benchmark:
        "../benchmarks/fastqc_multiqc-r2_benchmark.txt",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """
# ----- rule ----- #
