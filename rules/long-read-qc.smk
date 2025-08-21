logger.info('Running NanoPlot for quality control on long-read sequencing data')
rule NanoPlot:
    input:
        fastq = os.path.join(config["long_read_raw_data_path"], "{long_sample}.fq.gz"),
    output:
        report_html = "../01.qc/nanoplot/{long_sample}/NanoPlot-report.html",
    conda:
        "../envs/Nanoplot.yaml"
    log:
        "../logs/nanoplot/{long_sample}.nanoplot.log",
    params:
        report_dir = "../01.qc/nanoplot/{long_sample}",
        # nanoplot = config["software"]["qc"]["nanoplot"],
        title = "{long_sample}_nanoplot",
    message:
        "Running NanoPlot on {wildcards.long_sample}",
    threads: 
        config["threads"]["nanoplot"],
    shell:
        """
        NanoPlot -t {threads} \
                 --N50 \
                 --dpi 800 \
                 --fastq {input.fastq} \
                 --title {params.title} -o {params.report_dir} > {log} 2>&1        
        """

logger.info('Running fastplong for quality control on long-read sequencing data')
rule fastplong:
    input:
        fastq = os.path.join(config["long_read_raw_data_path"], "{long_sample}.fq.gz"),
    output:
        r1_trimmed = "../01.qc/long_read_trim/{long_sample}.fastplong.fq.gz",
        html_report = "../01.qc/long_read_trim/{long_sample}.fastplong.html",
        json_report = "../01.qc/long_read_trim/{long_sample}.fastplong.json",
    conda:
        "../envs/long-read-qc.yaml",
    log:
        "../logs/long_read_trim/{long_sample}.fastplong.log",
    params:
        quality_threshold = config["long_read_qc"]["qualified_quality_phred"],
        length_required = config["long_read_qc"]["length_required"],
    message:
        "Running fastplong on {input.fastq}",
    threads: 
        config["threads"]["fastplong"],
    shell:
        """
        fastplong -i {input.fastq} \
              -o {output.r1_trimmed} \
              --thread {threads} --length_required {params.length_required} \
              --qualified_quality_phred {params.quality_threshold} \
              -h {output.html_report}  -j {output.json_report} > {log} 2>&1
        """