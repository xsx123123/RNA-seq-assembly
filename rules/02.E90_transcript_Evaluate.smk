#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule E90_salmon_index:
    input:
        filtered = "../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa",
    output:
        salmon_dir = directory("../04.E90_transcript_Evaluate/salmon_index/transcripts_index"),
    conda:
        "../envs/salmon.yaml",
    log:
        "../logs/E90_transcript_Evaluate/salmon_index.log",
    benchmark:
        "../benchmarks/E90_transcript_Evaluate_salmon_index.txt",
    threads:
        config["threads"]["salmon_index"],
    shell:
        """
        salmon index -t {input.filtered} \
                     -i {output.salmon_dir} \
                     -p {threads} &> {log}
        """

rule E90_salmon_mapping:
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
    benchmark:
        "../benchmarks/E90_transcript_Evaluate_{sample}_salmon_quant.txt",
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

rule E90_abundance_estimates_to_matrix:
    input:
        salmon_quant = expand("../04.E90_transcript_Evaluate/salmon_quant/{sample}/quant.sf", sample=load_samples.keys()),
    output:
        transrate_salmon = "../04.E90_transcript_Evaluate/abundance_estimates_to_matrix/salmon_quant.isoform.counts.matrix",
    conda:
        "../envs/Trinity.yaml",
    log:
        "../logs/E90_transcript_Evaluate/abundance_estimates_to_matrix.log",
    benchmark:
        "../benchmarks/abundance_estimates_to_matrix.txt",
    params:
        transrate_salmon_prefix = '../04.E90_transcript_Evaluate/abundance_estimates_to_matrix/salmon_quant',
        method = config['abundance_estimates']['est_method'],
    threads:1
    shell:
        """
        abundance_estimates_to_matrix.pl --est_method {params.method} \
                                         --gene_trans_map  none \
                                         --name_sample_by_basedir \
                                         --out_prefix {params.transrate_salmon_prefix} \
                                         {input.salmon_quant} &> {log}
        """

rule E90_ExN50_filtered_length_dedup:
    input:
        transrate_salmon = '../04.E90_transcript_Evaluate/abundance_estimates_to_matrix/salmon_quant.isoform.counts.matrix',
        transcript = "../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa",
    output:
        ExN50_result = "../04.E90_transcript_Evaluate/ExN50_filtered/ExN50.transcript.stats",
    conda:
        "../envs/Trinity.yaml",
    log:
        "../logs/ExN50/ExN50_filtered_length_dedup.log",
    benchmark:
        "../benchmarks/ExN50_filtered_length_dedup.txt",
    threads:1
    shell:
        """
        contig_ExN50_statistic.pl {input.transrate_salmon} \
                                  {input.transcript} transcript | tee {output.ExN50_result} &> {log}
        """

rule E90_ExN50_filtered_length_dedup_visablity:
    input:
        ExN50_result = "../04.E90_transcript_Evaluate/ExN50_filtered/ExN50.transcript.stats",
    output:
        ExN50_plot = '../04.E90_transcript_Evaluate/ExN50_filtered/ExN50_transcript.pdf',
        ExN50_plot_png = '../04.E90_transcript_Evaluate/ExN50_filtered/ExN50_transcript.png',
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/ExN50/E90_ExN50_transcript_stats_plot.log",
    benchmark:
        "../benchmarks/E90_ExN50_transcript_stats_plot.txt",
    params:
        output_prefix = '../04.E90_transcript_Evaluate/ExN50_filtered/ExN50_transcript',
    threads:1
    shell:
        """
        python3 ./scripts/ExN50.py -i {input.ExN50_result} \
                         -o {params.output_prefix} \
                         -F png pdf &> {log}
        """
# ----- rule ----- #