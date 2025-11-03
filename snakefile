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
include: 'rules/00.common.smk'
include: 'rules/00.log.smk'
include: 'rules/00.id_convert.smk'
include: 'rules/00.get_all_input_dirs.smk'
include: 'rules/01.file_convert_md5.smk'
include: 'rules/01.short_read_qc.smk'
include: 'rules/01.Contamination_check.smk'
include: 'rules/01.short_read_clean.smk'
include: 'rules/01.long-read-qc_clean.smk'
include: 'rules/02.rnabloom_qc.smk'
include: 'rules/02.E90_filtered.smk'
include: 'rules/02.E90_transcript_Evaluate.smk'
include: 'rules/02.assembly_qc.smk'
include: 'rules/02.transcript2cds.smk'
include: 'rules/03.annotation_by_Diamond.smk'
include: 'rules/03.annotation_HMMER.smk'
include: 'rules/03.annotation_by_interproscan.smk'
include: 'rules/03.annotation_signal_peptides.smk'
# include: 'rules/04.RNA-seq-assembly-DEG.smk'
# --------- target rule --------- #
rule all:
    input:
        rna_assembly(config = config)
# --------- judge dependencies --------- #
judge_file_optimized(config = config)
# --------- target rule --------- #