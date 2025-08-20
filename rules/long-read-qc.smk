logger.info('Running NanoPlot for quality control on long-read sequencing data')
rule NanoPlot:
    input:
        fastq = os.path.join(config["long_read_raw_data_path"], "{long_sample}.fq.gz"),
    output:
        report_html = "../01.qc/nanoplot/{long_sample}/NanoPlot-report.html",
    conda:
        "../envs/qc.yaml"
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
        nanoplot -t {threads} \
                 --N50 \
                 --dpi 800 \
                 --fastq {input.fastq} \
                 --title {params.title} -o {params.report_dir} > {log} 2>&1        
        """ 