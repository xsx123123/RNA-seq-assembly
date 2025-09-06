rule rnabloom_long:
    input:
        long_reads = expand("../01.qc/long_read_trim/{sample}_long.fastq.gz", sample=long_read_samples.keys()),
        short_r1 = expand("../01.qc/trim/{sample}.R1.fastp.fq.gz", sample=load_samples.keys()),
        short_r2 = expand("../01.qc/trim/{sample}.R2.fastp.fq.gz", sample=load_samples.keys()),
    output:
        bloom = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.fa",
    conda:
        "../envs/rnabloom.yaml",
    params:
        outdir="../03.assembly/rnabloom_assembly",
    threads:
        config["threads"]["rnabloom"],
    shell:
        """
        rnabloom -long {input.long_reads} \
                 -t {threads} \
                 -outdir {params.outdir} \
                 -sef {input.short_r1} \
                 -sef {input.short_r2} \
                 -fpr 0.005 -overlap 200 -length 150 -lrop 0.7 -p 0.7 -lrrd 3
        """


rule transcripts_filter:
    input:
        bloom = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.fa",
    output:
        filtered = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.filtered.fa",
    conda:        
        "../envs/seqtk.yaml",                
    params:
        min_length = config["rnabloom"]["min_length"],
    log:
        "../logs/rnabloom/transcripts_filter.log",
    shell:
        """
        seqtk seq -L {params.min_length} \
              -A {input.bloom} > {output.filtered} 2> {log}
        """

rule transcripts_qc_busco:
    input:
        filtered = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.filtered.fa",
    output:
        busco_dir = directory("../03.assembly/rnabloom_assembly/busco/"),
    conda:
        "../envs/busco.yaml",
    log:
        "../logs/rnabloom/transcripts_qc_busco.log",
    params:
        lineage = config["busco"]["lineage"],
        mode = config["busco"]["mode"]"transcriptome",
    threads:
        config["threads"]["busco"],
    shell:
        """
        busco -i {input.filtered} \
              -o {output.busco_dir} \
              -l {params.lineage} \
              -m {params.mode} \
              -c {threads} 2> {log}
        """

rule transcripts_cluster:
    input:
        filtered = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.filtered.fa",
    output:
        clustered = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.filtered.dedup.fa",
    conda:
        "../envs/cd-hit.yaml",
    params:
        coverage = config["cdhit"]["coverage"],
    log:
        "../logs/rnabloom/transcripts_cluster.log",
    shell:
        """
        cd-hit-est -i {input.filtered} \
                   -o {output.clustered} \
                   -c {params.coverage} \
                   -n 11 -M 0 -T 22 -aL 0.9 -aS 0.9 -G 0
        """ 

rule transcripts_cluter_qc_busco:
    input:
        filtered = "../03.assembly/rnabloom_assembly/rnabloom.transcripts.filtered.dedup.fa",
    output:
        busco_dir = directory("../03.assembly/rnabloom_assembly/busco_clutser/"),
    conda:
        "../envs/busco.yaml",
    log:
        "../logs/rnabloom/transcripts_clutser_qc_busco.log",
    params:
        lineage = config["busco"]["lineage"],
        mode = config["busco"]["mode"]"transcriptome",
    threads:
        config["threads"]["busco"],
    shell:
        """
        busco -i {input.filtered} \
              -o {output.busco_dir} \
              -l {params.lineage} \
              -m {params.mode} \
              -c {threads} 2> {log}
        """