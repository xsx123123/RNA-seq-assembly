# Author : jian zhang 
# *---utf-8---*
# Version: 1.0v
# Author : jzhang
# ------- snakemake version check ------- #
from snakemake.utils import min_version
min_version("9.9.0")
# --------- main snakefile --------- #
configfile: "config/config.yaml"
configfile: "config.yaml"
# include all rules from the rules directory
include: 'rules/log.smk'
include: 'rules/id_convert.smk'
# include: 'rules/qc.smk'
# include: 'rules/trim.smk'
# include: 'rules/long-read-qc.smk'
# include: 'rules/STAR.smk'
# include: 'rules/rnabloom.smk'
# include: 'rules/Corset.smk'
include: 'rules/salmon.smk'
include: 'rules/transrate.smk'
include: 'rules/target.smk'