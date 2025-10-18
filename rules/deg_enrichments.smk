#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule DEG_by_salmon:
    input:
        filtered = "../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa",
    output:
        salmon_dir = directory("../04.E90_transcript_Evaluate/salmon_index/transcripts_index"),
    conda:
        "../envs/salmon.yaml",
    log:
        "../logs/E90_transcript_Evaluate/salmon_index.log",
    threads:
        config["threads"]["salmon_index"],
    shell:
        """
        salmon index -t {input.filtered} \
                     -i {output.salmon_dir} \
                     -p {threads} &> {log}
        """

rule Enrichments:
    input:
        r1 = "../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
        r2 = "../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz", 
        salmon_dir = "../04.E90_transcript_Evaluate/salmon_index/transcripts_index",
    output:
        quant = "../04.E90_transcript_Evaluate/salmon_quant/{sample}/quant.sf",
    conda:
        "../envs/salmon.yaml",
    log:
        "../logs/E90_transcript_Evaluate/{sample}_salmon_quant.log",
    threads:
        config["threads"]["salmon_quant"],
    params:
        quant = "../04.E90_transcript_Evaluate/salmon_quant/{sample}/",
        librarytype=config["salmon"]["librarytype"],
    shell:
        """
        salmon quant -i {input.salmon_dir} \
             -l {params.librarytype} \
             -1 {input.r1} \
             -2 {input.r2} \
             -p {threads} \
             --validateMappings \
             -o {params.quant} &> {log}
        """
# ----- rule ----- #