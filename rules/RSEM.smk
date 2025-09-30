rule rsem_quant:
    input:
        r1 = "../01.qc/trim/{sample}.R1.fastp.fq.gz",
        r2 = "../01.qc/trim/{sample}.R2.fastp.fq.gz", 
    output:
        isoforms="../04.rsem/{sample}.isoforms.results",
        genes="../04.rsem/{sample}.genes.results",
    params:
        ref=config["reference"],
        ref_prefix=config["rsem_bowtie2_index"],
        out_prefix="../04.rsem/{sample}" 
    conda:
        "../envs/rsem.yaml"
    threads: 6
    log:
        "../logs/rsem/{sample}.log"
    shell:
        """
        rsem-calculate-expression \
            --paired-end \
            --bowtie2 \
            --num-threads {threads} \
            --append-names \
            --estimate-rspd \
            {input.r1} {input.r2} \
            {params.ref_prefix} \
            {params.out_prefix} &> {log}
        """