# Set the default target to run
rule all:
    input:
        # Short-Read raw data QC and trimming reports
        # "../01.qc/multiqc/",
        # "../01.qc/multiqc_trim/",
        # expand("../01.qc/trim/{sample}.R1.fastp.fq.gz", sample=load_samples.keys()),
        # expand("../01.qc/trim/{sample}.R2.fastp.fq.gz", sample=load_samples.keys()),
        # expand("../01.qc/trim/{sample}.fastp.html", sample=load_samples.keys()),
        # expand("../01.qc/trim/{sample}.fastp.json", sample=load_samples.keys()),
        # Long-Read raw data QC and trißmming reports
        # expand("../01.qc/nanoplot/{long_sample}/NanoPlot-report.html",long_sample=long_read_samples.keys()),
        # expand("../01.qc/long_read_trim/{long_sample}.fastplong.fq.gz",long_sample=long_read_samples.keys()),
        # expand("../01.qc/long_read_trim/{long_sample}.fastplong.html",long_sample=long_read_samples.keys()),
        # expand("../01.qc/long_read_trim/{long_sample}.fastplong.json",long_sample=long_read_samples.keys()),
        # STAR mapping and BAM index
        # expand("../02.star/{sample}/{sample}.sort.bam",sample=load_samples.keys()),
        # expand("../02.star/{sample}/{sample}.sort.bam.bai",sample=load_samples.keys()),
        # "../02.star/multiqc_star/",
        # RNA-Bloom assembly
        # "../03.assembly/rnabloom_long_assembly"
        # RNA-Bloom transcripts redundancy removal
        expand("../04.mapping/salmon_quant/{sample}/",sample=load_samples.keys()),
        "../04.mapping/all_samples.TPM.matrix",
        "../04.mapping/E90_transcript_ids.txt",
        "../04.mapping/E90_transcript_busco/",
        "../04.mapping/E90_transcript_transrate/",