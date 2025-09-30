#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #

signalp -f short -n signalp.out Trinity.fasta.transdecoder.pep

tmhmm --short < Trinity.fasta.transdecoder.pep > tmhmm.out


rule Diamond_blast:
    input:
        pep = '../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa.TD2.pep',
    output:
        matches = '../05.transcript_annotation/TD2_pep_matches.tsv',
    conda:
        "../envs/td2.yaml",
    log:
        "../logs/diamond/diamond_blastp.log",
    threads:
        config["threads"]["TD2.LongOrfs"],
    params:
        swissprot = "/data/jzhang/reference/blast/swissprot/swissprot"
    shell:
        """
        diamond blastp -d {params.swissprot} \
                       -q {input.pep} \
                       -o {output.matches}  &> {log}
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
        python3 ../scripts/annotate_blast_results.py --blast_result {input.matches} \
                --uniprot_map {params.annotate_blast_uniport} \
                --output {output.annotated}  &> {log}
        """
# ----- rule ----- #