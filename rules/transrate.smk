rule transrate_assessment:
    input:
        transcript = "../04.mapping/rnabloom.transcripts.filtered.dedup_E90_transcript.fa",
        short_r1 = expand("../01.qc/trim/{sample}.R1.fastp.fq.gz", sample=load_samples.keys()),
        short_r2 = expand("../01.qc/trim/{sample}.R2.fastp.fq.gz", sample=load_samples.keys()),
    output:
        transrate_dir = directory("../04.mapping/E90_transcript_transrate/"),
    conda:
        "../envs/transrate.yaml",
    log:
        "../logs/salmon/E90_transcript_transrate.log",
    threads:
        config["threads"]["transrate"],
    shell:
        """
        transrate --assembly {input.transcript} \
                  --left "{",".join(input.short_r1)}" \
                  --right "{",".join(input.short_r2)}" \
                  --output {output.transrate_dir} \
                  --threads {threads} 2> {log}
        """