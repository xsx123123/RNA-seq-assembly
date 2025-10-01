#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule pep_signalp:
    input:
        pep = '../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa.TD2.pep',
    output:
        signalp_result = '../05.transcript_annotation/TD2_pep_signalp.out',
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/transcript_annotation/signalp.log",
    params:
        signalp = config['software']['signalp'],
    threads: 1
    shell:
        """
        {params.signalp} -f short \
                        -n {output.signalp_result} \
                        {input.input} &> {log}
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