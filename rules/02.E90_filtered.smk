#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
# ----- rule ----- #
rule E90_filtered_salmon_index:
    input:
        filtered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup.fa",
    output:
        salmon_dir = directory("../03.E90_filter/salmon_index/transcripts_index"),
    conda:
        "../envs/salmon.yaml",
    log:
        "../logs/salmon/E90_filtered_salmon_index.log",
    benchmark:
        "../benchmarks/E90_filtered_salmon_index.txt",
    threads:
        config["threads"]["salmon_index"],
    shell:
        """
        salmon index -t {input.filtered} \
                     -i {output.salmon_dir} \
                     -p {threads} &> {log}
        """

rule salmon_mapping:
    input:
        r1 = "../01.qc/short_read_trim/{sample}.R1.fastp.fq.gz",
        r2 = "../01.qc/short_read_trim/{sample}.R2.fastp.fq.gz", 
        salmon_dir = "../03.E90_filter/salmon_index/transcripts_index",
    output:
        quant = "../03.E90_filter/salmon_quant/{sample}/quant.sf",
    conda:
        "../envs/salmon.yaml",
    log:
        "../logs/salmon/{sample}_salmon_quant.log",
    benchmark:
        "../benchmarks/{sample}_salmon_quant.txt",
    threads:
        config["threads"]["salmon_quant"],
    params:
        quant = "../03.E90_filter/salmon_quant/{sample}/",
        librarytype=config["salmon"]["librarytype"],
    shell:
        """
        salmon quant -i {input.salmon_dir} \
             -l {params.librarytype} \
             -1 {input.r1} \
             -2 {input.r2} \
             -p {threads} \
             --validateMappings \
             -o {params.quant} &> {log}
        """

rule abundance_estimates_to_matrix:
    input:
        salmon_quant = expand("../03.E90_filter/salmon_quant/{sample}/quant.sf", sample=load_samples.keys()),
    output:
        transrate_salmon = "../03.E90_filter/abundance_estimates_to_matrix/salmon_quant.isoform.counts.matrix",
    conda:
        "../envs/Trinity.yaml",
    log:
        "../logs/abundance_estimates_to_matrix/abundance_estimates_to_matrix.log",
    benchmark:
        "../benchmarks/abundance_estimates_to_matrix.txt",
    params:
        transrate_salmon_prefix = '../03.E90_filter/abundance_estimates_to_matrix/salmon_quant',
        method = config['abundance_estimates']['est_method'],
    threads:1
    shell:
        """
        abundance_estimates_to_matrix.pl --est_method {params.method} \
                                         --gene_trans_map  none \
                                         --name_sample_by_basedir \
                                         --out_prefix {params.transrate_salmon_prefix} \
                                         {input.salmon_quant} &> {log}
        """

rule ExN50_filtered_length_dedup:
    input:
        transrate_salmon = "../03.E90_filter/abundance_estimates_to_matrix/salmon_quant.isoform.counts.matrix",
        transcript = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup.fa",
    output:
        ExN50_result = "../03.E90_filter/ExN50_filtered/ExN50.transcript.stats",
    conda:
        "../envs/Trinity.yaml",
    log:
        "../logs/ExN50/ExN50_filtered_length_dedup.log",
    benchmark:
        "../benchmarks/ExN50_filtered_length_dedup.txt",
    threads:1
    shell:
        """
        contig_ExN50_statistic.pl {input.transrate_salmon} \
                                  {input.transcript} transcript | tee {output.ExN50_result} &> {log}
        """

rule ExN50_filtered_length_dedup_visablity:
    input:
        ExN50_result = "../03.E90_filter/ExN50_filtered/ExN50.transcript.stats",
    output:
        ExN50_plot = '../03.E90_filter/ExN50_filtered/ExN50_transcript.pdf',
        ExN50_plot_png = '../03.E90_filter/ExN50_filtered/ExN50_transcript.png',
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/ExN50/ExN50_transcript_stats_plot.log",
    benchmark:
        "../benchmarks/ExN50_transcript_stats_plot.txt",
    params:
        output_prefix = '../03.E90_filter/ExN50_filtered/ExN50_transcript',
    threads:1
    shell:
        """
        python3 ./scripts/ExN50.py -i {input.ExN50_result} \
                         -o {params.output_prefix} \
                         -F png pdf &> {log}
        """

rule E90_aggregate_quants:
    input:
        quants = expand("../03.E90_filter/salmon_quant/{sample}/quant.sf", sample=load_samples.keys()),
    output:
        tpm_matrix = "../03.E90_filter/all_samples.TPM.matrix"
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/salmon/aggregate_salmon_tpm.log"
    benchmark:
        "../benchmarks/aggregate_salmon_tpm.txt",
    shell:
        """
        python3 scripts/merge_salmon_quants.py {input.quants} \
                                               -o {output.tpm_matrix} &> {log}
        """

rule E90_core_transcripts:
    input:
        quant = "../03.E90_filter/all_samples.TPM.matrix",
    output:
        counts = "../03.E90_filter/E90_transcript_ids.txt",
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/salmon/E90_filter.log",
    benchmark:
        "../benchmarks/E90_filter.txt",
    threads: 1
    shell:
        """
        python3 ./scripts/E90_filter.py {input.quant} \
                                        0.9 {output.counts} &> {log}
        """

rule E90_Extert_core_transcripts:
    input:
        E90_transcript = "../03.E90_filter/E90_transcript_ids.txt",
        filtered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup.fa",
    output:
        transcript = "../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa",
    conda:
        "../envs/seqtk.yaml",
    log:
        "../logs/salmon/seqtk_E90_filter.log",
    benchmark:
        "../benchmarks/seqtk_E90_filter.txt",
    threads: 1
    shell:
        """
        seqtk subseq {input.filtered} \
                     {input.E90_transcript} > {output.transcript} 2> {log}
        """

rule E90_Exterted_transcripts:
    input:
        E90_transcript_ids = "../03.E90_filter/E90_transcript_ids.txt",
        filtered = "../02.assembly/rnabloom_assembly/rnabloom.transcripts.length_filtered.dedup.fa",
    output:
        remove_transcript = "../03.E90_filter/remove_transcript.fa",
    conda:
        "../envs/seqtk.yaml",
    log:
        "../logs/salmon/remove_transcript.log",
    benchmark:
        "../benchmarks/remove_transcript.txt",
    threads: 1
    shell:
        """
        awk 'NR==FNR {{ids[$1]=1; next}} \
             /^>/ {{id=substr($1,2); if (id in ids) {{p=0}} else {{p=1}}}} p' \
             {input.E90_transcript_ids} {input.filtered} > {output.remove_transcript} 2> {log}
        """

rule busco_assessment_for_E90:
    input:
        transcript = "../03.E90_filter/rnabloom.transcripts.length_filtered.dedup_E90_transcript.fa",
    output:
        busco_dir = directory("../03.E90_filter/length_filtered_cd-hit-est_E90_transcript_busco/"),
    conda:
        "../envs/busco.yaml",
    log:
        "../logs/salmon/E90_transcript_busc.log",
    benchmark:
        "../benchmarks/E90_transcript_busc.txt",
    params:
        lineage = config["busco"]["lineage"],
        mode = config["busco"]["mode"],
    threads:
        config["threads"]["busco_salmon"],
    shell:
        """
        busco -i {input.transcript} \
              -o {output.busco_dir} \
              --lineage_dataset {params.lineage} \
              -m {params.mode} \
              -c {threads} &> {log}
        """

rule busco_assessment_for_remove:
    input:
        transcript = "../03.E90_filter/remove_transcript.fa",
    output:
        busco_dir = directory("../03.E90_filter/length_filtered_cd-hit-est_remove_transcrip_busco/"),
    conda:
        "../envs/busco.yaml",
    log:
        "../logs/salmon/remove_transcrip_busco.log",
    benchmark:
        "../benchmarks/remove_transcrip_busco.txt",
    params:
        lineage = config["busco"]["lineage"],
        mode = config["busco"]["mode"],
    threads:
        config["threads"]["busco_salmon"],
    shell:
        """
        busco -i {input.transcript} \
              -o {output.busco_dir} \
              --lineage_dataset {params.lineage} \
              -m {params.mode} \
              -c {threads} &> {log}
        """
# ----- rule ----- #