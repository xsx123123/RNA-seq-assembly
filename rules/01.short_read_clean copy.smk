#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
rule short_read_fastp:
    input:
        md5_check = "../01.qc/md5_check.tsv",
    output:
        r1_trimmed = "../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
        r2_trimmed = "../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz",
        html_report = "../01.qc/short_read_trim/{sample}.fastp.html",
        json_report = "../01.qc/short_read_trim/{sample}.fastp.json",
    conda:
        "../envs/fastp.yaml",
    log:
        "../logs/01.short_read_trim/{sample}.fastp.log",
    message:
        "Running Fastp on {wildcards.sample} r1 and {wildcards.sample} r2",
    benchmark:
        "../benchmarks/{sample}_fastp_benchmark.txt",
    params:
        length_required = config["trim"]["length_required"],
        quality_threshold = config["trim"]["quality_threshold"],
        adapter_fasta = config["trim"]["adapter_fasta"],
        r1 = os.path.join(config["raw_data_path"],
                          config['convert_md5'],
                          "{sample}",
                          "{sample}" + config['r1_suffix']),
        r2 = os.path.join(config["raw_data_path"],
                          config['convert_md5'],
                          "{sample}",
                          "{sample}" + config['r2_suffix']),
    threads: 
        config["threads"]["fastp"],
    shell:
        """
        fastp -i {params.r1} -I {params.r2} \
              -o {output.r1_trimmed} -O {output.r2_trimmed} \
              --thread {threads} \
              --length_required  {params.length_required} \
              --qualified_quality_phred {params.quality_threshold} \
              --adapter_fasta {params.adapter_fasta} \
              -g -V \
              -h {output.html_report} \
              -j {output.json_report} &> {log}
        """

rule multiqc_trim:
    input:
        md5_check = "../01.qc/md5_check.tsv",
        fastp_report = expand("../01.qc/short_read_trim/{sample}.fastp.html", sample=samples.keys()),
    output:
        report_dir = directory("../01.qc/multiqc_short_read_trim/")
    conda:
        "../envs/multiqc.yaml",
    message:
        "Running MultiQC to aggregate fastp reports",
    benchmark:
        "../benchmarks/multiqc_fastp_benchmark.txt",
    params:
        fastqc_reports = "../01.qc/short_read_trim/",
        report = "multiqc_short_read_trim_report.html",
        title = "short_read_trim-multiqc-report",
    log:
        "../logs/01.short_read_trim/multiqc_trim.log",
    shell:
        """
        multiqc {params.fastqc_reports} \
                --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} &> {log}
        """
# ----- rule ----- #