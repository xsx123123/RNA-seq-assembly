#!/usr/bin/env python3
from rich import print as rich_print
# Set the default target to run
def rna_assembly(config:dict = None) -> list:
    # short-read + long-read hybrid RNA-seq assembly
    hybrid_rna_assembly = [
            "../01.qc/short_read_r1_multiqc/",
            "../01.qc/short_read_r2_multiqc/"
        ]
    if config['long_read']:
        hybrid_rna_assembly.extend(expand("../01.qc/long_read_qc/{long_sample}/NanoPlot-report.html",
                                          long_sample=long_read_samples.keys()))
        hybrid_rna_assembly.extend(expand("../01.qc/long_read_trim/{long_sample}.fastplong.fq.gz",
                                          long_sample=long_read_samples.keys()))
        hybrid_rna_assembly.extend(expand("../01.qc/long_read_trim_qc/{long_sample}/NanoPlot-report.html",
                                          long_sample=long_read_samples.keys()))
        #hybrid_rna_assembly.append("../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup/") 
        #hybrid_rna_assembly.extend(expand("../03.E90_filter/salmon_quant/{sample}/quant.sf",
        #                                    sample=load_samples.keys()))

    rich_print(hybrid_rna_assembly)

    return  hybrid_rna_assembly
