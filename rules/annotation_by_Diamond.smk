#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule Diamond_blast:
    input:
        pep = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.pep',
    output:
        matches = '../05.transcript_annotation/TD2_pep_matches.tsv',
    conda:
        "../envs/td2.yaml",
    log:
        "../logs/diamond/diamond_blastp.log",
    threads:
        config["threads"]["diamond_blastp_threads"],
    params:
        swissprot = "/data/jzhang/reference/blast/swissprot/swissprot"
    shell:
        """
        diamond blastp --threads {threads} \
                       -d {params.swissprot} \
                       -q {input.pep} \
                       -o {output.matches} 2> {log}
        """

rule uniport_ann:
    input:
        matches = '../05.transcript_annotation/TD2_pep_matches.tsv',
    output:
        annotated = '../05.transcript_annotation/TD2_pep_matches_annotated.tsv',
    params:
        annotate_blast_uniport = config["annotate"]["uniport_database"],
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/transcript_annotation/annotate_blast_results.log",
    shell:
        """
        python3 ./scripts/annotate_blast_results.py --blast_result {input.matches} \
                --uniprot_map {params.annotate_blast_uniport} \
                --output {output.annotated} 2> {log}
        """

rule uniport_go:
    input:
        matches_annotated = '../05.transcript_annotation/TD2_pep_matches_annotated.tsv',
    output:
        go_dataset = '../05.transcript_annotation/TD2_pep_matches_annotated_go.tsv',
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/transcript_annotation/annotate_blast_results_go_dataset.log",
    shell:
        """
        python3 ./scripts/extract_go_terms.py \
                --input {input.matches_annotated} \
                --output {output.go_dataset} 2> {log}
        """
# ----- rule ----- #