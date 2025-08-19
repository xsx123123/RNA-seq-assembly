rule qc:
    input:
        r1 = os.path.join(config["raw_data_path"], "{sample}.R1.fq.gz"),
        r2 = os.path.join(config["raw_data_path"], "{sample}.R2.fq.gz"),
    output:
        r1_html = "../01.qc/fastqc/{sample}.R1_fastqc.html",
        r1_zip = "../01.qc/fastqc/{sample}.R1_fastqc.zip", 
        r2_html = "../01.qc/fastqc/{sample}.R2_fastqc.html",
        r2_zip = "../01.qc/fastqc/{sample}.R2_fastqc.zip", 
    log:
        r1 = "../logs/qc/{sample}.r1.fastqc.log",
        r2 = "../logs/qc/{sample}.r2.fastqc.log",
    params:
        fastqc = config["software"]["qc"]["fastqc"],
        out_dir = "../01.qc/fastqc",
    message:
        "Running FastQC on {input.r1} and {input.r2}",
    threads: 1
    shell:
        """
        {params.fastqc} {input.r1} -o {params.out_dir} --threads {threads} > {log.r1} 2>&1
        {params.fastqc} {input.r2} -o {params.out_dir} --threads {threads} > {log.r2} 2>&1
        """

rule multiqc:
    input:
        fastqc_reports = expand("../01.qc/fastqc/{sample}.R1_fastqc.zip", sample=load_samples.keys()) + \
                         expand("../01.qc/fastqc/{sample}.R2_fastqc.zip", sample=load_samples.keys())
    output:
        report = "../01.qc/multiqc/multiqc_report.html",
    message:
        "Running MultiQC to aggregate FastQC reports",
    params:
        out_dir = "../01.qc/multiqc",
    log:
        "../logs/multiqc.log",
    shell:
        """
        multiqc {input.fastqc_reports} -o {params.out_dir} -n multiqc_report.html > {log} 2>&1
        """