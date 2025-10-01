#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule interproscan:
    input:
        pep = '../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa.TD2.pep',
    output:
        clean_pep = '../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa.clean.pep',
        ann = directory('../05.transcript_annotation/TD2_pep_interproscan_annotation/'),
    log:
        "../logs/transcript_annotation/interproscan.log",
    threads:
        config["threads"]["interproscan"],
    params:
        interproscan = config['software']['interproscan']
    shell:
        """
        sed sed 's/\*//g' {input.pep} > {output.clean_pep} &&
        {params.interproscan} -i {output.clean_pep} \
                              -f tsv \
                              -d {output.ann} \
                              -cpu {threads} \
                              --goterms \
                              &> {log}
        """
# ----- rule ----- #