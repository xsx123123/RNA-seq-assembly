#!/usr/bin/env python3
# *---utf-8---*
# Version: 0.1.2v
# Author : jzhang
# ------- snakemake version check ------- #
from snakemake.utils import min_version
min_version("9.9.0")
# --------- main snakefile --------- #
configfile: "config/config.yaml"
configfile: "config.yaml"
# --------- snakemake rule --------- #
# include all rules from the rules directory
include: 'rules/common.smk'
include: 'rules/log.smk'
include: 'rules/id_convert.smk'
include: 'rules/short_read_qc.smk'
include: 'rules/short_read_clean.smk'
include: 'rules/long-read-qc_clean.smk'
include: 'rules/rnabloom_qc.smk'
include: 'rules/E90_filtered.smk'
include: 'rules/E90_transcript_Evaluate.smk'
include: 'rules/assembly_qc.smk'
include: 'rules/transcript2cds.smk'
include: 'rules/annotation_by_Diamond.smk'
include: 'rules/annotation_HMMER.smk'
include: 'rules/RNA-seq-assembly-DEG.smk'
# --------- target rule --------- #
rule all:
    input:
        rna_assembly(config = config)
# --------- target rule --------- #