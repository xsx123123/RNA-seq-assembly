rule fastp:
    input:
        r1 = os.path.join(config["raw_data_path"], "{sample}.R1.fq.gz"),
        r2 = os.path.join(config["raw_data_path"], "{sample}.R2.fq.gz"),
    output:
        r1_trimmed = "../01.qc/trim/{sample}.R1.fastp.fq.gz",
        r2_trimmed = "../01.qc/trim/{sample}.R2.fastp.fq.gz",
        html_report = "../01.qc/trim/{sample}.fastp.html",
        json_report = "../01.qc/trim/{sample}.fastp.json",
    log:
        "../logs/trim/{sample}.fastp.log",
    message:
        "Running Fastp on {input.r1} and {input.r2}",
    params:
        fastp = config["software"]["qc"]["fastp"],
        length_required = config["trim"]["length_required"],
        quality_threshold   = config["trim"]["quality_threshold"],
    threads: 4
    shell:
        """
        {params.fastp} -i {input.r1} -I {input.r2}  
        -o {output.r1_trimmed} -O {output.r2_trimmed} 
        --thread {threads} --length_required  {params.length_required}
        --qualified_quality_phred {params.quality_threshold} -g -V
        -h {output.html_report}  -j {output.json_report} > {log} 2>&1
        """