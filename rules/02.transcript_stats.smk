#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule transcript_stats:
    input:
        stats_bloom = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.fa",
        stats_filtered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.fa",
        stats_clustered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup.fa",
        stats_transcript = "../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa",
        stats_remove_transcript = "../03.E90_filter/remove_transcript.fa",
        stats_cds = '../05.transcript_annotation/rnabloom_transcript_LongOrfs.fa.TD2.cds',
    output:
        stats_transcript = "../03.E90_filter/transcript_stats.tsv",
    conda:
        "../envs/seqkit.yaml",
    log:
        "../logs/03.E90_filter/seqkit_stats.log",
    benchmark:
        "../benchmarks/seqkit_stats.txt",
    threads:
        1
    shell:
        """
        seqkit stats {input.stats_bloom} \
                     {input.stats_filtered} \
                     {input.stats_clustered} \
                     {input.stats_transcript} \
                     {input.stats_remove_transcript} \
                     {input.stats_cds} > {output.stats_transcript} 2>{log}
        """
# ----- rule ----- #