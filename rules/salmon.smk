rule salmon_index:
    input:
        filtered = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.filtered.dedup.fa",
    output:
        salmon_dir = directory("../04.mapping/salmon_index/transcripts_index"),
    conda:
        "../envs/salmon.yaml",
    log:
        "../logs/salmon/salmon_index.log",
    threads:
        config["threads"]["salmon_index"],
    shell:
        """
        salmon index -t {input.filtered} \
                     -i {output.salmon_dir} \
                     -p {threads} 2> {log}
        """

rule salmon_mapping:
    input:
        r1 = "../01.qc/trim/{sample}.R1.fastp.fq.gz",
        r2 = "../01.qc/trim/{sample}.R2.fastp.fq.gz", 
        salmon_dir = directory("../04.mapping/salmon_index/transcripts_index"),
    output:
        quant = directory("../04.mapping/salmon_quant/{sample}/"),
    conda:
        "../envs/salmon.yaml",
    log:
        "../logs/salmon/{sample}_salmon_quant.log",
    threads:
        config["threads"]["salmon_quant"],
    params:
        librarytype=config["salmon"]["librarytype"],
    shell:
        """
        salmon quant -i {input.salmon_dir} \
             -l {params.librarytype} \
             -1 {input.r1} \
             -2 {input.r2} \
             -p {threads} \
             --validateMappings \
             -o {output.quant} 2> {log}
        """

rule aggregate_quants:
    input:
        quants = expand("../04.mapping/salmon_quant/{sample}/quant.sf", sample=load_samples.keys())
    output:
        tpm_matrix = "../04.mapping/all_samples.TPM.matrix"
    params:

    shell:
        """
        ../scripts/aggregate_salmon_tpm.py {input.quants} {output.tpm_matrix}
        """

rule core_transcripts:
    input:
        quant = "../04.mapping/all_samples.TPM.matrix",
    output:
        counts = "../04.mapping/E90_transcript_ids.txt",
    conda:
        "../envs/corset.yaml",
    log:
        "../logs/corset/corset_counts.log",
    threads: 1
    shell:
        """
        ../scripts/E90_filter.py {input.quant} 0.9  {output.counts}
        """