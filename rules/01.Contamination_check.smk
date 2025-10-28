#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
rule short_read_fastq_screen_r1:
    input:
        md5_check = "../01.qc/md5_check.tsv",
    output:
        fastq_screen_result = "../01.qc/fastq_screen_r1/{sample}_R1_screen.txt",
    log:
        "../logs/01.short_read_qc_r1/{sample}.r1.fastq_screen.log",
    params:
        out_dir = "../01.qc/fastq_screen_r1/",
        fastq_screen_dir = config['fastq_screen']['path'],
        conf = config['fastq_screen']['conf'],
        subset = config['fastq_screen']['subset'],
        aligner = config['fastq_screen']['aligner'],
        r1 = os.path.join(config["raw_data_path"],
                          config['convert_md5'],
                          "{sample}",
                          "{sample}" + config['r1_suffix']),
    message:
        "Running fastq_screen on {wildcards.sample} r1",
    benchmark:
        "../benchmarks/{sample}_r1_fastq_screen_benchmark.txt",
    threads: 
        config['threads']['fastq_screen'],
    shell:
        """
        {params.fastq_screen_dir} --threads  {threads} \
                     --force \
                     --subset  {params.subset} \
                     --aligner  {params.aligner} \
                     --conf {params.conf} \
                     --outdir {params.out_dir} \
                     {params.r1} &> {log}
        """

rule short_read_fastq_screen_r2:
    input:
        md5_check = "../01.qc/md5_check.tsv",
    output:
        fastq_screen_result = "../01.qc/fastq_screen_r2/{sample}_R2_screen.txt",
    log:
        "../logs/01.short_read_qc_r2/{sample}.r2.fastq_screen.log",
    params:
        out_dir = "../01.qc/fastq_screen_r2/",
        conf = config['fastq_screen']['conf'],
        subset = config['fastq_screen']['subset'],
        aligner = config['fastq_screen']['aligner'],
        r2 = os.path.join(config["raw_data_path"],
                          config['convert_md5'],
                          "{sample}",
                          "{sample}" + config['r2_suffix']),
    message:
        "Running fastq_screen on {wildcards.sample} r2",
    benchmark:
        "../benchmarks/{sample}_r2_fastq_screen_benchmark.txt",
    threads: 
        config['threads']['fastq_screen'],
    shell:
        """
        fastq_screen --threads  {threads} \
                     --force \
                     --subset  {params.subset} \
                     --aligner  {params.aligner} \
                     --conf {params.conf} \
                     --outdir {params.out_dir} \
                     {params.r2} &> {log}
        """

rule fastq_screen_multiqc_r1:
    input:
        fastqc_files_r1 = expand("../01.qc/fastq_screen_r1/{sample}_R1_screen.txt", sample=load_samples.keys()),
    output:
        report_dir = directory("../01.qc/fastq_screen_multiqc_r1/")
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to aggregate R1 fastq screen reports",
    params:
        fastqc_reports = "../01.qc/fastq_screen_r1",
        report = "multiqc_r1_fastq_screen_report.html",
        title = "r1-fastq-screen-multiqc-report",
    log:
        "../logs/01.multiqc/multiqc-fastq-screen-r1.log",
    benchmark:
        "../benchmarks/fastqc_multiqc-fastq-screen-r1_benchmark.txt",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """

rule fastq_screen_multiqc_r2:
    input:
        fastqc_files_r1 = expand("../01.qc/fastq_screen_r2/{sample}_R2_screen.txt", sample=load_samples.keys()),
    output:
        report_dir = directory("../01.qc/fastq_screen_multiqc_r2/")
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to aggregate R2 fastq screen reports",
    params:
        fastqc_reports = "../01.qc/fastq_screen_r2",
        report = "multiqc_r2_fastq_screen_report.html",
        title = "r2-fastq-screen-multiqc-report",
    log:
        "../logs/01.multiqc/multiqc-fastq-screen-r2.log",
    benchmark:
        "../benchmarks/fastqc_multiqc-fastq-screen-r2_benchmark.txt",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """
# ----- rule ----- #
