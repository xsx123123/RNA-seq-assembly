#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule RNA_seq_assembly_DEG:
    input:
        salmon_quant = "../05.assembly/salmon_quant",
        sample_info = "./data/sample_info.csv",
        transcript_ids = "./data/E90_transcript_ids.txt",
    output:
        deg_result = "../06.DEG/deg_result.csv",
    conda:
        "../envs/rna-seq_R.yaml",
    log:
        "../logs/DEG/DEG-analysis.log",
    params:
        pval = config['deg']['pval'],
        lfc = config['deg']['lfc'],
    shell:
        """
        Rscript ./scripts/run_deg_analysis.R \
        --salmon_dir {input.salmon_quant} \
        --sample_info {input.sample_info} \
        --transcript_ids {input.transcript_ids} \
        --output_dir  {output.deg_result} \
        --pval {params.pval} \
        --lfc {params.lfc} &> {log}
        """
# ----- rule ----- #