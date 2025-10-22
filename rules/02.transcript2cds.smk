#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule transcript2cds:
    input:
        transcript_fasta = "../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa",
    output:
        cds_fasta = "../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript/longest_orfs.cds",
        pep_fasta = "../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript/longest_orfs.pep",
        gff = "../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript/longest_orfs.gff3", 
    conda:
        "../envs/td2.yaml",
    log:
        "../logs/transcript2cds/transcript2cds.log",
    benchmark:
        "../benchmarks/transcript2cds.txt",
    threads:
        config["threads"]["TD2.LongOrfs"],
    params:
        link_dir = "../05.transcript_annotation/",
        link_data = "../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa",
        TD2_LongOrfs_dir = '../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript/',
    shell:
        """
        cp {input.transcript_fasta} {params.link_dir} &&
        TD2.LongOrfs -t {params.link_data} \
                     --threads {threads} \
                     -O {params.TD2_LongOrfs_dir} &>{log}
        """

rule pep_search_by_mmseqs:
    input:
        pep_fasta = "../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript/longest_orfs.pep",
    output:
        mmseqs_result = "../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript/longest_orfs_alnRes.m8",
    conda:
        "../envs/mmseqs2.yaml",
    log:
        "../logs/pep_search_by_mmseq/pep_search_by_mmseq.log",
    benchmark:
        "../benchmarks/pep_search_by_mmseq.txt",
    params:
        mmseqs_database = config["mmseqs"]["swissprot_database"],
        temp_dir = "../05.transcript_annotation/tmp",
    threads:
        config["threads"]["mmseqs"],
    shell:
        """
        mmseqs easy-search {input.pep_fasta} \
                           {params.mmseqs_database} \
                           {output.mmseqs_result} \
                           {params.temp_dir} \
                           -s 7.0 \
                           --threads {threads} &>{log}
        """

rule TD2Predict_Integrating_homology:
    input:
        mmseqs_result = "../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript/longest_orfs_alnRes.m8",
    output:
        pep = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.pep',
        cds = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.cds',
        gff3 = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.gff3',
        bed = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.bed',
    conda:
        "../envs/td2.yaml",
    log:
        "../logs/transcript2cds/transcript2cds.log",
    benchmark:
        "../benchmarks/transcript2cds.txt",
    threads:
        config["threads"]["TD2.LongOrfs"],
    params:
        TD2_LongOrfs_dir = "../05.transcript_annotation/rnabloom.transcripts.length_filtered.dedup_E90_transcript",
        TD2_convert_dir = "../05.transcript_annotation/rnabloom_transcript",
        TD2_analysis_dir = "../05.transcript_annotation/",
        TD2_transcript = './rnabloom_transcript',
    shell:
        """
        cp -r {params.TD2_LongOrfs_dir} {params.TD2_convert_dir} &&
        cp -r {params.TD2_LongOrfs_dir}.fa {params.TD2_convert_dir}.fa &&
        cd {params.TD2_analysis_dir} &&
        TD2.Predict -t {params.TD2_transcript}.fa \
                    --retain-mmseqs-hits {input.mmseqs_result} -v &> {log}
        """

rule TD2_busco:
    input:
        transcript = '../05.transcript_annotation/rnabloom_transcript.fa.TD2.cds',
    output:
        busco_dir = directory("../05.transcript_annotation/TD2_CDS_busco/"),
    conda:
        "../envs/busco.yaml",
    log:
        "../logs/transcript2cds/TD2_CDS_busco.log",
    benchmark:
        "../benchmarks/TD2_CDS_busco.txt",
    params:
        lineage = config["busco"]["lineage"],
        mode = config["busco"]["mode"],
    threads:
        config["threads"]["busco_salmon"],
    shell:
        """
        busco -i {input.transcript} \
              -o {output.busco_dir} \
              --lineage_dataset {params.lineage} \
              -m {params.mode} \
              -c {threads} &> {log}
        """
# ----- rule ----- #