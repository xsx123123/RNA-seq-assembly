#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule transrate_assessment:
    input:
        short_r1 = expand("../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz", sample=load_samples.keys()),
        short_r2 = expand("../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz", sample=load_samples.keys()),
    output:
        transrate_dir = directory("../05.transcript_annotation/E90_transcript_transrate/"),
    conda:
        "../envs/transrate.yaml",
    log:
        "../logs/transcript_annotation/E90_transcript_transrate.log",
    benchmark:
        "../benchmarks/E90_transcript_transrate.txt",
    params:
        short_r1 = ",".join(expand("../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz", sample=load_samples.keys())),
        short_r2 = ",".join(expand("../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz", sample=load_samples.keys())),
        transcript = "../04.mapping/rnabloom.transcripts.filtered.dedup_E90_transcript.fa",
    threads:
        config["threads"]["transrate"],
    shell:
        """
        transrate --assembly {params.transcript} \
                  --left {params.short_r1} \
                  --right {params.short_r2} \
                  --output {output.transrate_dir} \
                  --threads {threads} &> {log}
        """
# ----- rule ----- #