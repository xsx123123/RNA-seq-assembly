#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
# logger.info('Running NanoPlot for quality control on long-read sequencing data')
rule NanoPlot:
    input:
        md5_check = '../01.qc/md5_check_short_read.tsv',
        fastq = os.path.join('../00.raw_data',
                              config['short_read_convert_md5'],
                              "{long_sample}/{long_sample}.fq.gz"),
    output:
        report_html = "../01.qc/long_read_qc/{long_sample}/NanoPlot-report.html",
    conda:
        "../envs/Nanoplot.yaml"
    log:
        "../logs/nanoplot/{long_sample}.nanoplot.log",
    benchmark:
        "../benchmarks/{long_sample}_NanoPlot-report.txt",
    params:
        report_dir = "../01.qc/long_read_qc/{long_sample}",
        # nanoplot = config["software"]["qc"]["nanoplot"],
        title = "{long_sample}_nanoplot",
    message:
        "Running NanoPlot on {wildcards.long_sample}",
    threads: 
        config["threads"]["nanoplot"],
    shell:
        """
        NanoPlot -t {threads} \
                 --N50 \
                 --dpi 800 \
                 --fastq {input.fastq} \
                 --title {params.title} \
                 -o {params.report_dir} &> {log}        
        """

# logger.info('Running fastplong for quality control on long-read sequencing data')
rule fastplong_trim:
    input:
        md5_check = '../01.qc/md5_check_short_read.tsv',
        fastq = os.path.join('../00.raw_data',
                              config['short_read_convert_md5'],
                              "{long_sample}/{long_sample}.fq.gz"),
    output:
        r1_trimmed = "../01.qc/long_read_trim/{long_sample}.fastplong.fq.gz",
        html_report = "../01.qc/long_read_trim/{long_sample}.fastplong.html",
        json_report = "../01.qc/long_read_trim/{long_sample}.fastplong.json",
    conda:
        "../envs/long-read-qc.yaml",
    log:
        "../logs/long_read_trim/{long_sample}.fastplong.log",
    benchmark:
        "../benchmarks/{long_sample}_fastplong.txt",
    params:
        quality_threshold = config["long_read_qc"]["qualified_quality_phred"],
        length_required = config["long_read_qc"]["length_required"],
    message:
        "Running fastplong on {input.fastq}",
    threads: 
        config["threads"]["fastplong"],
    shell:
        """
        fastplong -i {input.fastq} \
                  -o {output.r1_trimmed} \
                  --thread {threads} --length_required {params.length_required} \
                  --qualified_quality_phred {params.quality_threshold} \
                  -h {output.html_report}  -j {output.json_report} &> {log}
        """

# logger.info('Running NanoPlot for quality control on trim long-read sequencing data')
rule NanoPlot_trim:
    input:
        fastq = "../01.qc/long_read_trim/{long_sample}.fastplong.fq.gz",
    output:
        report_html = "../01.qc/long_read_trim_qc/{long_sample}/NanoPlot-report.html",
    conda:
        "../envs/Nanoplot.yaml"
    log:
        "../logs/nanoplot/{long_sample}.trim.nanoplot.log",
    benchmark:
        "../benchmarks/{long_sample}_trim.nanoplot.txt",
    params:
        report_dir = "../01.qc/long_read_trim_qc/{long_sample}",
        title = "{long_sample}_trim_nanoplot",
    message:
        "Running NanoPlot on Trim {wildcards.long_sample}",
    threads: 
        config["threads"]["nanoplot"],
    shell:
        """
        NanoPlot -t {threads} \
                 --N50 \
                 --dpi 800 \
                 --fastq {input.fastq} \
                 --title {params.title} \
                 -o {params.report_dir} &> {log}        
        """
# ----- rule ----- #