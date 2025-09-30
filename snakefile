#!/usr/bin/env python3
# *---utf-8---*
# Version: 1.0v
# Author : jzhang
# ------- snakemake version check ------- #
from snakemake.utils import min_version
min_version("9.9.0")
# --------- main snakefile --------- #
configfile: "config/config.yaml"
configfile: "config.yaml"
# --------- snakemake rule --------- #
# include all rules from the rules directory
include: 'rules/target.smk'
include: 'rules/log.smk'
include: 'rules/id_convert.smk'
include: 'rules/short_read_qc.smk'
include: 'rules/long-read-qc_clean.smk'
include: 'rules/rnabloom_qc.smk'
include: 'rules/salmon_E90.smk'
# --------- target rule --------- #
rule all:
    input:
        rna_assembly(config = config)
# --------- target rule --------- #