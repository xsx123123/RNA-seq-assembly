#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule Trinity:
    input:
        short_r1 = expand("../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
                            sample=load_samples.keys()),
        short_r2 = expand("../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz",
                            sample=load_samples.keys()),
    output:
        trinity_fasta = "../02.trinity/{sample}.Trinity.fasta",
        trinity_transcripts = "../02.trinity/{sample}.Trinity.transcripts.fasta",
        trinity_log = "../logs/trinity/{sample}.Trinity.log",
    conda:
        "../envs/trinity.yaml",
    log:
        "../logs/trinity/{sample}.Trinity.log",
    params:
        min_kmer_cov = config["trinity"]["min_kmer_cov"],
        max_memory = config["trinity"]["max_memory"],
    message:
        "Running Trinity on {input.r1} and {input.r2}",
    threads: 
        config["threads"]["trinity"],
    shell:
        """
        trinity --seqType fq \
                --left {input.short_r1} \
                --right {input.short_r2} \
                --max_memory {params.max_memory} \
                --min_kmer_cov {params.min_kmer_cov} \
                --output {output.trinity_fasta} \
                --CPU {threads} \
                --log {output.trinity_log} &> {log}
        """
# ----- rule ----- #
