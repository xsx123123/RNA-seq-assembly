#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule pep_signalp:
    input:
        pep = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.pep',
    output:
        signalp_result = directory('../05.transcript_annotation/TD2_pep_signalp'),
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/transcript_annotation/signalp.log",
    params:
        signalp = config['software']['signalp'],
        run_mode = config['signalp']['mode'],
    threads: config['threads']['signalp'],
    shell:
        """
        {params.signalp} --fastafile {input.pep} \
                         --organism other \
                         --output_dir {output.signalp_result} \
                         --format txt \
                         --mode {params.run_mode} \
                          &> {log}
        """

rule pep_tmhmm:
    input:
        pep = '../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa.TD2.pep',
    output:
        tmhmm_result = '../05.transcript_annotation/TD2_pep_tmhmm.out',
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/transcript_annotation/tmhmm.log",
    threads: 1
    params:
        tmhmm = config['software']['tmhmm'],
    shell:
        """
        {params.tmhmm} --short < {input.pep} > {output.tmhmm_result} &> {log}
        """
# ----- rule ----- #