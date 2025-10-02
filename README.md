*Author*  : jzhang \n
*Email*   : zhangjian199567@outlook.com \n
*Date*    : 2025-10-2 \n
*Version* : 0.1.2v
## Introduce
### `RNA-seq-assembly`
`RNA-seq-assembly`是基于`short-read` + `ONT long-read`的转录组组装流程，包含数据质控、转录本组装、转录本注释、差异分析与富集分析的`snakmake`流程。
### Workflow
`RNA-seq-assembly`转录本组装流程主要由以下步骤组成:
raw-data qc -> raw-data clean -> hybrid assembly by rnabloom -> assembly qc -> Transcript Redundancy Reduction -> Transcript Annotation -> DEG -> GO/KEGG enrichment
## Dependencies & Install
`snakemake = 9.9.0`
```bash
pip3 install pandas
pip3 install rich
pip3 install loguru
```
## Run
### Install `miniconda3`
```bash 
# Install miniconda3
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```
### Install `snakemake=9.9.0`
```bash
# fix conda channel_priority -> flexible
conda config --set channel_priority flexible
# Create snakemake conda environment and install snakemak
conda create -n snakemake
# activate snakemake
conda activate snakemake
# install snakemake=9.9.0
conda install bioconda::snakemake=9.9.0
```
### Run de novo transcriptome assembly by short-read & long-read
```bash
# Creating conda environments via mamba
snakemake --use-conda --conda-create-envs-only --conda-frontend mamba
# Clean install conda/mamba envs
snakemake --use-conda --conda-cleanup-envs
# Runing the snakemake pipeline via mamba
snakemake --cores=30 -p --conda-frontend mamba --use-conda
```
## Reference