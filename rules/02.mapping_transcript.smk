#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule bowtie_build:
    input:
        transcript_fasta = "../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa",
    output:
        bowtie_index_bt2_1 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.1.bt2",
        bowtie_index_bt2_2 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.2.bt2",
        bowtie_index_bt2_3 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.3.bt2",
        bowtie_index_bt2_4 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.4.bt2",
        bowtie_index_bt2_rev_1 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.rev.1.bt2",
        bowtie_index_bt2_rev_2 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.rev.2.bt2",
    conda:
        "../envs/bowtie2.yaml",
    log:
        "../logs/bowtie_index/bowtie_index.log",
    benchmark:
        "../benchmarks/bowtie_index.txt",
    threads:
        config["threads"]["bowtie2_build"],
    params:
        bowtie2_build_index_dir = '../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index',
    shell:
        """
        bowtie2-build  --threads {threads} \
                       {input.transcript_fasta} \
                       {params.bowtie2_build_index_dir} &>{log}
        """

rule bowtie_mapping:
    input:
        bowtie_index_bt2_1 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.1.bt2",
        bowtie_index_bt2_2 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.2.bt2",
        bowtie_index_bt2_3 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.3.bt2",
        bowtie_index_bt2_4 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.4.bt2",
        bowtie_index_bt2_rev_1 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.rev.1.bt2",
        bowtie_index_bt2_rev_2 = "../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index.rev.2.bt2",
        r1 = "../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
        r2 = "../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz", 
    output:
        sam = temp('../04.E90_transcript_Evaluate/bowtie2_mapping/{sample}.sam'),
    conda:
        "../envs/bowtie2.yaml",
    log:
        "../logs/bowtie2_mapping/bowtie2_{sample}.log",
    benchmark:
        "../benchmarks/bowtie2_{sample}.txt",
    params:
        bowtie2_build_index_dir = '../04.E90_transcript_Evaluate/bowtie2_build/bowtie2_build_index',
    threads:
        config["threads"]["bowtie_mapping"],
    shell:
        """
        bowtie2 -x {params.bowtie2_build_index_dir} \
                -1 {input.r1} \
                -2 {input.r2}  \
                -S {output.sam} -p {threads} &>{log}
        """

rule convert_sort_index:
    input:
        sam = '../04.E90_transcript_Evaluate/bowtie2_mapping/{sample}.sam',
    output:
        bam = '../04.E90_transcript_Evaluate/bowtie2_mapping/{sample}.bam',
        sort_bam = '../02.mapping/bwa_mem2/sort_index/{sample}.sort.bam',
        sort_bam_bai = '../02.mapping/bwa_mem2/sort_index/{sample}.sort.bam.bai',
    conda:
        "../envs/bwa2.yaml",
    message:
        "Running samtools sort & index for {wildcards.sample}",
    log:
        "../logs/02.mapping/bwa_sort_index_{sample}.log",
    benchmark:
            "../benchmarks/{sample}_bam_sort_index_benchmark.txt",
    threads: 
        config["threads"]["samtools_threads"],
    shell:
        """
        (samtools samtools view -@ {threads} -Sbh -o {output.bam} {input.sam} &&
        samtools sort -@ {threads} -o {output.sort_bam} {output.bam} &&
        samtools index -@ {threads} {output.sort_bam}) 2>{log}
        """

# ----- rule ----- #