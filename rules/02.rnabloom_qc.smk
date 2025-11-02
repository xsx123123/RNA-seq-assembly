#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ----- rule ----- #
rule rnabloom_long:
    input:
        long_reads = expand("../01.qc/long_read_trim/{sample}.fastplong.fq.gz",
                            sample=long_read_samples.keys()),
        short_r1 = expand("../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
                            sample=load_samples.keys()),
        short_r2 = expand("../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz",
                            sample=load_samples.keys()),
    output:
        bloom = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.fa",
    conda:
        "../envs/rnabloom.yaml",
    log:
        "../logs/rnabloom/rnabloom_long.log",
    benchmark:
        "../benchmarks/rnabloom_long.txt",
    params:
        outdir="../02.assembly/rnabloom_assembly",
    threads:
        config["threads"]["rnabloom"],
    shell:
        """
        rnabloom -long {input.long_reads} \
                 -t {threads} \
                 -outdir {params.outdir} \
                 -sef {input.short_r1} \
                 -sef {input.short_r2} &> {log}
        """

rule transcripts_filter:
    input:
        bloom = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.fa",
    output:
        filtered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.fa",
    conda:        
        "../envs/seqtk.yaml",                
    params:
        min_length = config["rnabloom"]["min_length"],
    log:
        "../logs/rnabloom/transcripts_filter.log",
    benchmark:
        "../benchmarks/transcripts_filter.txt",
    shell:
        """
        seqtk seq -L {params.min_length} -A {input.bloom} > {output.filtered} 2> {log}
        """

rule transcripts_qc_busco:
    input:
        filtered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.fa",
    output:
        busco_dir = directory("../02.assembly/rnabloom_assembly/length_filtered_busco/"),
    conda:
        "../envs/busco.yaml",
    log:
        "../logs/rnabloom/transcripts_length_filtered_busco.log",
    benchmark:
        "../benchmarks/transcripts_length_filtered_busco.txt",
    params:
        lineage = config["busco"]["lineage"],
        mode = config["busco"]["mode"]
    threads:
        config["threads"]["busco"],
    shell:
        """
        busco -i {input.filtered} \
              -o {output.busco_dir} \
              -l {params.lineage} \
              -m {params.mode} \
              -c {threads} &> {log}
        """

rule transcripts_cluster:
    input:
        filtered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.fa",
    output:
        clustered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup.fa",
    conda:
        "../envs/cd-hit.yaml",
    params:
        coverage = config["cdhit"]["coverage"],
    log:
        "../logs/rnabloom/transcripts_dedup_cd-hit.log",
    benchmark:
        "../benchmarks/transcripts_dedup_cd-hit.txt",
    shell:
        """
        cd-hit-est -i {input.filtered} \
                   -o {output.clustered} \
                   -c {params.coverage} \
                   -n 11 -M 0 -T 22 -aL 0.9 -aS 0.9 -G 0 &> {log}
        """ 

rule transcripts_cluter_qc_busco:
    input:
        filtered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup.fa",
    output:
        busco_dir = directory("../02.assembly/rnabloom_assembly/length_filtered_cd-hit-est_busco/"),
    conda:
        "../envs/busco.yaml",
    log:
        "../logs/rnabloom/rnabloom.transcripts.length_filtered.dedup.log",
    benchmark:
        "../benchmarks/rnabloom.transcripts.length_filtered.dedup.txt",
    params:
        lineage = config["busco"]["lineage"],
        mode = config["busco"]["mode"],
    threads:
        config["threads"]["busco"],
    shell:
        """
        busco -i {input.filtered} \
              -o {output.busco_dir} \
              -l {params.lineage} \
              -m {params.mode} \
              -c {threads} &> {log}
        """
# ----- rule ----- #