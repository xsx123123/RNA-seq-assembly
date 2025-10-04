#!/usr/bin/env python3
from rich import print as rich_print
# Set the default target to run
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
    # long-read qc & trim & clean result
    if config['long_read']:
        hybrid_rna_assembly.extend(expand("../01.qc/long_read_qc/{long_sample}/NanoPlot-report.html",
                                          long_sample=long_read_samples.keys()))
        hybrid_rna_assembly.extend(expand("../01.qc/long_read_trim/{long_sample}.fastplong.fq.gz",
                                          long_sample=long_read_samples.keys()))
        hybrid_rna_assembly.extend(expand("../01.qc/long_read_trim_qc/{long_sample}/NanoPlot-report.html",
                                          long_sample=long_read_samples.keys()))
    # RNA-Assembly result
    # RNA-Assembly busco result
    hybrid_rna_assembly.append("../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup/") 
    hybrid_rna_assembly.extend(expand("../03.E90_filter/salmon_quant/{sample}/quant.sf",
                                        sample=load_samples.keys()))
    # Fliter E90_transcript
    hybrid_rna_assembly.append('../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa')
    # Fliter E90_transcript busco result
    hybrid_rna_assembly.append("../03.E90_filter/E90_transcript_busco/")
    # Evaluate ExN50 for E90_transcrip
    hybrid_rna_assembly.extend(expand("../04.E90_transcript_Evaluate/salmon_quant/{sample}/quant.sf",
                                          sample=load_samples.keys()))
    # TD2.Predict CDS
    # hybrid_rna_assembly.append("../05.transcript_annotation/E90_transcript_transrate/")
    hybrid_rna_assembly.append('../05.transcript_annotation/rnabloom_transcript.fa.TD2.pep')
    hybrid_rna_assembly.append('../05.transcript_annotation/rnabloom_transcript.fa.TD2.cds')
    hybrid_rna_assembly.append('../05.transcript_annotation/rnabloom_transcript.fa.TD2.gff3')
    hybrid_rna_assembly.append('../05.transcript_annotation/rnabloom_transcript.fa.TD2.bed')
    # Trinotate annotation
    # --- diamond swissprot annotation --- #
    hybrid_rna_assembly.append("../05.transcript_annotation/TD2_pep_matches_annotated.tsv")
    # --- PFAM --- #
    hybrid_rna_assembly.append('../05.transcript_annotation/TrinotatePFAM.out')
    # interproscan annotation
    hybrid_rna_assembly.append("../05.transcript_annotation/rnabloom_transcript.fa.TD2.pep.clean")
    hybrid_rna_assembly.append("../05.transcript_annotation/TD2_pep_interproscan_annotation/")
    # Print Target rule           
    if config['print_target']:
        rich_print(hybrid_rna_assembly)
    return  hybrid_rna_assembly
