logger.info('Running Trinity for de novo transcriptome assembly')



rule Trinity:
    input:
        r1 = os.path.join(config["raw_data_path"], "{sample}.R1.fq.gz"),
        r2 = os.path.join(config["raw_data_path"], "{sample}.R2.fq.gz"),
    output:
        trinity_fasta = "../02.trinity/{sample}.Trinity.fasta",
        trinity_transcripts = "../02.trinity/{sample}.Trinity.transcripts.fasta",
        trinity_log = "../logs/trinity/{sample}.Trinity.log",
    conda:
        "../envs/trinity.yaml",
    log:
        "../logs/trinity/{sample}.Trinity.log",
    params:
        trinity = config["software"]["trinity"],
        min_kmer_cov = config["trinity"]["min_kmer_cov"],
        max_memory = config["trinity"]["max_memory"],
    message:
        "Running Trinity on {input.r1} and {input.r2}",
    threads: 
        config["threads"]["trinity"],
    shell:
        """
        {params.trinity} --seqType fq \
                         --left {input.r1} \
                         --right {input.r2} \
                         --max_memory {params.max_memory} \
                         --min_kmer_cov {params.min_kmer_cov} \
                         --output {output.trinity_fasta} \
                         --CPU {threads} \
                         --log {output.trinity_log} > {log} 2>&1
        """ 
