logger.info('Running STAR for RNA-seq read mapping to reference genome')

rule STAR_mapping:
     input:
          r1 = "../01.qc/trim/{sample}.R1.fastp.fq.gz",
          r2 = "../01.qc/trim/{sample}.R2.fastp.fq.gz",
     output:
          bam = "../02.star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
          sj = "../02.star/{sample}/{sample}.SJ.out.tab",
     conda:
          "../envs/mapping.yaml",
     log:
          "../logs/star/{sample}.STAR.log",
     params:
          star_index = config["star_index"],
          out_prefix = "../02.star/{sample}/{sample}.",
     message:
          "Running STAR mapping for {wildcards.sample}",
     threads: 
          config["threads"]["star_mapping"],
     shell:
          """
          STAR --runThreadN {threads} \
                           --genomeDir {params.star_index} \
                           --readFilesIn {input.r1} {input.r2} \
                           --readFilesCommand zcat \
                           --outFileNamePrefix {params.out_prefix} \
                           --outSAMtype BAM SortedByCoordinate \
                           --outSAMunmapped None \
                           --twopassMode Basic \
                           --limitBAMsortRAM 60000000000 &> {log}
          """

logger.info('Running samtools to sort and index BAM files')

rule index_bam:
     input:
          bam = "../02.star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
     output:
          sort = "../02.star/{sample}/{sample}.sort.bam",
          bai = "../02.star/{sample}/{sample}.sort.bam.bai",
     conda:
          "../envs/mapping.yaml",
     log:
          sort = "../logs/star/{sample}.index_bam.log",
          bai = "../logs/star/{sample}.index_bam.log",
     message:
          "Indexing BAM file for {wildcards.sample}",
     threads: 
          config["threads"]["samtools_sort"],
     shell:
          """
          samtools sort -@ {threads} {input.bam} -o {output.sort} &> {log.sort}
          samtools index -@ {threads} {output.sort} -o {output.bai} &> {log.bai}
          """

logger.info('Running multiqc to summarize STAR mapping logs')

rule multiqc_star:
     input:
          star_result = expand("../02.star/{sample}/{sample}.sort.bam", sample=load_samples.keys()),
     output:
          report_dir = directory("../02.star/multiqc_star/")
     conda:
          "../envs/qc.yaml",
     message:
          "Running MultiQC to aggregate STAR mapping logs",
     params:
          star_dir = "../02.star/",
          report_dir = "../02.star/multiqc_star/",
          report = "multiqc_star_report.html",
          title = "star-mapping-multiqc-report",
     log:
          "../logs/star/multiqc_star.log",
     shell:
          """
          multiqc {params.star_dir} \
                  --outdir {output.report_dir} \
                  -i {params.title} \
                  -n {params.report} &> {log}
          """