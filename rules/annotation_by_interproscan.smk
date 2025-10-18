#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule interproscan:
    input:
        pep = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.pep',
    output:
        clean_pep = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.pep.clean',
        interproscan_result = '../05.transcript_annotation/TD2_pep_interproscan_annotation/rnabloom_transcript.fa.TD2.pep.clean.tsv',
    log:
        "../logs/transcript_annotation/interproscan.log",
    threads:
        config["threads"]["interproscan"],
    params:
        interproscan = config['software']['interproscan'],
        ann_dir = '../05.transcript_annotation/TD2_pep_interproscan_annotation/',
    shell:
        r"""
        mkdir -p {params.ann_dir} &&
        sed 's/\*//g' {input.pep} > {output.clean_pep} &&
        {params.interproscan} -i {output.clean_pep} \
                              -f tsv \
                              -d {params.ann_dir} \
                              -cpu {threads} \
                              --goterms \
                              &> {log}
        """

rule interproscan_format:
    input:
        ann = '../05.transcript_annotation/TD2_pep_interproscan_annotation/rnabloom_transcript.fa.TD2.pep.clean.tsv',
    output:
        summary = "../05.transcript_annotation/TD2_pep_interproscan_annotation/interproscan_ann.summary",
        go = "../05.transcript_annotation/TD2_pep_interproscan_annotation/interproscan_ann.gopathway",
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/transcript_annotation/preprocess_interpro.log",
    threads:
        config['threads']['preprocess_interpro'],
    shell:
        """
        python3 ./scripts/preprocess_interpro.py \
                --summary_out {output.summary} \
                --go_map_out  {output.go} \
                {input.ann} \
                --cores {threads} &> {log}
        """
# ----- rule ----- #