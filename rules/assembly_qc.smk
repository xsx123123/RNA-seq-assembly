logger.info('Running transrate for quality control on rnabloom assessment data')
rule transrate_assessment:
    input:
        short_r1 = expand("../01.qc/trim/{sample}.R1.fastp.fq.gz", sample=load_samples.keys()),
        short_r2 = expand("../01.qc/trim/{sample}.R2.fastp.fq.gz", sample=load_samples.keys()),
    output:
        transrate_dir = directory("../05.transcript/E90_transcript_transrate/"),
    conda:
        "../envs/transrate.yaml",
    log:
        "../logs/salmon/E90_transcript_transrate.log",
    params:
        short_r1 = ",".join(expand("../01.qc/trim/{sample}.R1.fastp.fq.gz", sample=load_samples.keys())),
        short_r2 = ",".join(expand("../01.qc/trim/{sample}.R2.fastp.fq.gz", sample=load_samples.keys())),
        transcript = "../04.mapping/rnabloom.transcripts.filtered.dedup_E90_transcript.fa",
    threads:
        config["threads"]["transrate"],
    shell:
        """
        transrate --assembly {params.transcript} \
                  --left {params.short_r1} \
                  --right {params.short_r2} \
                  --output {output.transrate_dir} \
                  --threads {threads} &> {log}
        """

logger.info('Running TrinityStats for quality control on rnabloom assessment data')
rule TrinityStats_assembly:
    input:
        transcript = "../04.mapping/rnabloom.transcripts.filtered.dedup_E90_transcript.fa",
    output:
        transrate_dir = directory("../05.transcript/TrinityStats/"),
    conda:
        "../envs/Trinity.yaml",
    log:
        "../05.transcript/TrinityStats/TrinityStats.log",
    threads:1
    shell:
        """
        TrinityStats.pl {input.transcript} &> {log}
        """

logger.info('Running abundance_estimates_to_matrix for quality control on rnabloom assessment data')
rule abundance_estimates_to_matrix_filter:
    input:
        salmon_quant = expand("../04.mapping/salmon_quant/{sample}/quant.sf", sample=load_samples.keys()),
    output:
        transrate_salmon = "../05.transcript/abundance_estimates_to_matrix/salmon_quant.isoform.counts.matrix",
    conda:
        "../envs/Trinity.yaml",
    log:
        "../logs/abundance_estimates_to_matrix/abundance_estimates_to_matrix.log",
    params:
        transrate_salmon_prefix = '../05.transcript/abundance_estimates_to_matrix/salmon_quant',
        workdir = '../04.mapping/salmon_quant',
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


logger.info('Running ExN50 for quality control on rnabloom assessment data')
rule ExN50_filter:
    input:
        transrate_salmon = "../05.transcript/abundance_estimates_to_matrix/salmon_quant.isoform.counts.matrix",
        transcript = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.filtered.dedup.fa",
    output:
        ExN50_result = directory("../05.transcript/ExN50_filtered/")
    conda:
        "../envs/Trinity.yaml",
    log:
        "../logs/ExN50/ExN50.log",
    threads:1
    shell:
        """
        contig_ExN50_statistic.pl {input.transrate_salmon} \
                                  {input.transcript} &> {log}
        """
        