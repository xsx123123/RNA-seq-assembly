logger.info('Running FastQC for quality control on raw sequencing data')

rule qc:
    input:
        r1 = os.path.join(config["raw_data_path"], "{sample}.R1.fq.gz"),
        r2 = os.path.join(config["raw_data_path"], "{sample}.R2.fq.gz"),
    output:
        r1_html = "../01.qc/fastqc/{sample}.R1_fastqc.html",
        r1_zip = "../01.qc/fastqc/{sample}.R1_fastqc.zip", 
        r2_html = "../01.qc/fastqc/{sample}.R2_fastqc.html",
        r2_zip = "../01.qc/fastqc/{sample}.R2_fastqc.zip", 
    conda:
        "../envs/qc.yaml",
    log:
        r1 = "../logs/qc/{sample}.r1.fastqc.log",
        r2 = "../logs/qc/{sample}.r2.fastqc.log",
    params:
        # fastqc = config["software"]["qc"]["fastqc"],
        out_dir = "../01.qc/fastqc",
    message:
        "Running FastQC on {input.r1} and {input.r2}",
    threads: 1
    shell:
        """
        fastqc {input.r1} -o {params.out_dir} --threads {threads} > {log.r1} 2>&1 &&
        fastqc {input.r2} -o {params.out_dir} --threads {threads} > {log.r2} 2>&1
        """

logger.info('Run MultiQC to summarize fastqc QC reports')
rule multiqc:
    input:
        fastqc_files_r1 = expand("../01.qc/fastqc/{sample}.R1_fastqc.zip", sample=load_samples.keys()),
        fastqc_files_r2 = expand("../01.qc/fastqc/{sample}.R2_fastqc.zip", sample=load_samples.keys()),
    output:
        report_dir = directory("../01.qc/multiqc/")
    conda:
        "../envs/qc.yaml",
    message:
        "Running MultiQC to aggregate FastQC reports",
    params:
        fastqc_reports = "../01.qc/fastqc/",
        report = "multiqc_trim_report.html",
        # multiqc = config["software"]["qc"]["multiqc"],
        title = "raw-data-multiqc-report",
    log:
        "../logs/qc/multiqc.log",
    shell:
        """
        multiqc {params.fastqc_reports} --outdir {output.report_dir} \
                -i {params.title} \
                -n {params.report} > {log} 2>&1
        """