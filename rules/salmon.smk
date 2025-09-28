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
        salmon_dir = "../04.mapping/salmon_index/transcripts_index",
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
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/salmon/aggregate_salmon_tpm.log"
    shell:
        """
        python3 scripts/merge_salmon_quants.py {input.quants} -o {output.tpm_matrix} 2> {log}
        """

rule core_transcripts:
    input:
        quant = "../04.mapping/all_samples.TPM.matrix",
    output:
        counts = "../04.mapping/E90_transcript_ids.txt",
    conda:
        "../envs/python3.yaml",
    log:
        "../logs/salmon/E90_filter.log",
    threads: 1
    shell:
        """
        python3 ./scripts/E90_filter.py {input.quant} 0.9  {output.counts} 2> {log}
        """

rule Extert_core_transcripts:
    input:
        E90_transcript = "../04.mapping/E90_transcript_ids.txt",
        filtered = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.filtered.dedup.fa",
    output:
        transcript = "../04.mapping/rnabloom.transcripts.filtered.dedup_E90_transcript.fa",
    conda:
        "../envs/seqtk.yaml",
    log:
        "../logs/salmon/seqtk_E90_filter.log",
    threads: 1
    shell:
        """
        seqtk subseq {input.filtered} {input.E90_transcript} > {output.transcript}  2> {log}
        """

rule busco_assessment:
    input:
        transcript = "../04.mapping/rnabloom.transcripts.filtered.dedup_E90_transcript.fa",
    output:
        busco_dir = directory("../04.mapping/E90_transcript_busco/"),
    conda:
        "../envs/busco.yaml",
    log:
        "../logs/salmon/E90_transcript_busc.log",
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
              -c {threads} 2> {log}
        """