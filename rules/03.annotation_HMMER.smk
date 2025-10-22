#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule hmmscan:
    input:
        pep = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.pep',
    output:
        matches = '../05.transcript_annotation/TrinotatePFAM.out',
    conda:
        "../envs/hmmscan.yaml",
    log:
        "../logs/hmmscan/hmmscan_pep.log",
    benchmark:
        "../benchmarks/hmmscan_pep.txt",
    params:
        Pfam_A = config['annotate']['Pfam_A'],
    threads:
        config["threads"]["hmmscan"],
    shell:
        """
        hmmscan --cpu {threads} \
                --domtblout {output.matches} \
                {params.Pfam_A} \
                {input.pep} &> {log}
        """
# ----- rule ----- #