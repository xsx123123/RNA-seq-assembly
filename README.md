*Author*  : zhang jian
*Date*    : 2025-9-30
*Version* : 1.2v
## Introduce
### `RNA-seq-assembly`
`RNA-seq-assembly`是基于`short-read` + `ONT long-read`的转录组组装流程，包含数据质控、转录本组装、转录本注释、差异分析与富集分析的`snakmake`流程。
## Workflow


## Install
```bash
pip3 install pandas
pip3 install rich
pip3 install loguru
```

## Run
setup snakemake analysis pipeline by mamba
```bash
# fix conda channel_priority -> flexible
conda config --set channel_priority flexible
# Creating conda environments via mamba
snakemake --use-conda --conda-create-envs-only --conda-frontend mamba
# Clean install conda/mamba envs
snakemake --use-conda --conda-cleanup-envs
# Runing the snakemake pipeline via mamba
snakemake --cores=30 -p --conda-frontend mamba --use-conda
```