#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# loading packages
from loguru import logger
from pathlib import Path
from typing import Dict, Union
from rich import print as rich_print
# Target rule function
def rna_assembly(config:dict = None) -> list:
    """
    short-read & long-read RNA-seq assembly & annotation pipeline target
    """
    # short-read raw-data qc result
    hybrid_rna_assembly = [
            "../01.qc/short_read_r1_multiqc/",
            "../01.qc/short_read_r2_multiqc/"            
        ]
    # short-read trim & clean result
    hybrid_rna_assembly.extend(expand("../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
                                          sample=load_samples.keys()))
    hybrid_rna_assembly.extend(expand("../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz",
                                          sample=load_samples.keys()))
    hybrid_rna_assembly.extend(expand("../01.qc/short_read_trim/{sample}.fastp.html",
                                          sample=load_samples.keys()))
    hybrid_rna_assembly.extend(expand("../01.qc/short_read_trim/{sample}.fastp.json",
                                          sample=load_samples.keys()))
    hybrid_rna_assembly.append("../01.qc/multiqc_short_read_trim/")
    if config['fastq_screen']['run']:
        logger.info("fastq_screen analysis is enabled")
        hybrid_rna_assembly.append("../01.qc/fastq_screen_multiqc_r1/")
        hybrid_rna_assembly.append("../01.qc/fastq_screen_multiqc_r2/")
    else:
        logger.info("skipping fastq_screen analysis")
    # long-read qc & trim & clean result
    if config['long_read']:
        logger.info("long-read qc & trim & clean is enabled")
        hybrid_rna_assembly.extend(expand("../01.qc/long_read_qc/{long_sample}/NanoPlot-report.html",
                                          long_sample=long_read_samples.keys()))
        hybrid_rna_assembly.extend(expand("../01.qc/long_read_trim/{long_sample}.fastplong.fq.gz",
                                          long_sample=long_read_samples.keys()))
        hybrid_rna_assembly.extend(expand("../01.qc/long_read_trim_qc/{long_sample}/NanoPlot-report.html",
                                          long_sample=long_read_samples.keys()))
    else:
        logger.info("skipping long-read qc & trim & clean analysis")
    # RNA-Assembly result
    # RNA-Assembly busco result
    hybrid_rna_assembly.append("../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup/") 
    hybrid_rna_assembly.extend(expand("../03.E90_filter/salmon_quant/{sample}/quant.sf",
                                        sample=load_samples.keys()))
    # Fliter E90_transcript
    hybrid_rna_assembly.append('../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa')
    # Fliter E90_transcript busco result
    hybrid_rna_assembly.append("../03.E90_filter/E90_transcript_busco/")
    hybrid_rna_assembly.append("../03.E90_filter/remove_transcript.fa")
    hybrid_rna_assembly.append("../03.E90_filter/remove_transcrip_busco/")
    # Evaluate ExN50 for E90_transcrip
    hybrid_rna_assembly.extend(expand("../04.E90_transcript_Evaluate/salmon_quant/{sample}/quant.sf",
                                          sample=load_samples.keys()))
    hybrid_rna_assembly.append("../04.E90_transcript_Evaluate/abundance_estimates_to_matrix/salmon_quant.isoform.counts.matrix")
    hybrid_rna_assembly.append("../04.E90_transcript_Evaluate/ExN50_filtered/ExN50.transcript.stats")
    # TD2.Predict CDS
    # hybrid_rna_assembly.append("../05.transcript_annotation/E90_transcript_transrate/")
    hybrid_rna_assembly.append('../05.transcript_annotation/rnabloom_transcript_LongOrfs.fa.TD2.pep')
    hybrid_rna_assembly.append('../05.transcript_annotation/rnabloom_transcript_LongOrfs.fa.TD2.cds')
    hybrid_rna_assembly.append('../05.transcript_annotation/rnabloom_transcript_LongOrfs.fa.TD2.gff3')
    hybrid_rna_assembly.append('../05.transcript_annotation/rnabloom_transcript_LongOrfs.fa.TD2.bed')
    hybrid_rna_assembly.append("../05.transcript_annotation/TD2_CDS_busco/")
    # Trinotate annotation
    # --- diamond swissprot annotation --- #
    hybrid_rna_assembly.append("../05.transcript_annotation/TD2_pep_matches_annotated.tsv")
    # --- PFAM --- #
    hybrid_rna_assembly.append('../05.transcript_annotation/TrinotatePFAM.out')
    # interproscan annotation
    hybrid_rna_assembly.append("../05.transcript_annotation/rnabloom_transcript_LongOrfs.fa.TD2.clean")
    hybrid_rna_assembly.append("../05.transcript_annotation/TD2_pep_interproscan_annotation/rnabloom_transcript.fa.TD2.pep.clean.tsv")
    hybrid_rna_assembly.append("../05.transcript_annotation/TD2_pep_interproscan_annotation/interproscan_ann.summary")
    hybrid_rna_assembly.append("../05.transcript_annotation/TD2_pep_interproscan_annotation/interproscan_ann.gopathway")
    hybrid_rna_assembly.append("../05.transcript_annotation/TD2_pep_matches_annotated_go.tsv")
    # Prediction Signal peptide by SignalP 6.0
    # hybrid_rna_assembly.append('../05.transcript_annotation/TD2_pep_signalp')
    # Print Target rule           
    if config['print_target']:
        rich_print(hybrid_rna_assembly)
    return  hybrid_rna_assembly


def _check_file_existence(path_str: str, name: str) -> None:
    """
    help function to check file existence
    """
    path = Path(path_str)
    if not path.exists():
        error_msg = f'Place check the {name} path: {path_str}'
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)

def judge_file_optimized(config: Dict = None) -> None:
    """
    judge if the file is exist for config yaml file
    """
    if config is None:
        return

    software_paths = [
        (config["software"].get("tmhmm"), "tmhmm"),
        (config["software"].get("signalp"), "signalp"),
        (config["software"].get("interproscan"), "interproscan"),
    ]

    database_paths = [
        (config["busco"].get("lineage"), "BUSCO lineage database"),
        (config["mmseqs"].get("swissprot_database"), "swissprot database"),
        (config["annotate"].get("uniport_database"), "uniport database"),
        (config["annotate"].get("Pfam_A"), "Pfam A database"),
    ]

    all_paths = software_paths + database_paths

    for path_str, name in all_paths:
        if path_str:
            _check_file_existence(path_str, name)
