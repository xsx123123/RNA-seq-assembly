#!/usr/bin/snakemake
# -*- coding: utf-8 -*-
import os
# ------- rule ------- #
rule seq_preprocessor_short_read:
    input:
        md5 = config['raw_data_path'],
        # md5 = get_all_input_dirs(load_samples.keys(),config = config,raw_data_path_config = 'raw_data_path'),
    output:
        md5_check = directory(os.path.join('../00.raw_data',config['short_read_convert_md5'])),
        md5_check_json = os.path.join('../00.raw_data',config['short_read_convert_md5'],"raw_data_md5.json"),
        link_r1 = expand(os.path.join('../00.raw_data',
                                      config['short_read_convert_md5'],
                                      "{sample}/{sample}_R1.fq.gz"),
                                      sample=samples.keys()),
        link_r2 = expand(os.path.join('../00.raw_data',
                                      config['short_read_convert_md5'],
                                      "{sample}/{sample}_R2.fq.gz"),
                                      sample=samples.keys()),
    message:
        "Running seq_preprocessor on raw data data",
    benchmark:
        "../benchmarks/seq_preprocessor.txt",
    params:
        raw_data_path = config['raw_data_path'],
        md5 = config['short_read_md5'],
    log:
        "../logs/01.qc/seq_preprocessor.txt",
    threads: 1
    shell:
        """
        ./scripts/seq_preprocessor/target/release/seq_preprocessor -i  {params.raw_data_path} \
                -o {output.md5_check} \
                --library-type short-read \
                --md5-name {params.md5} \
                --json-report {output.md5_check_json} &> {log}
        """

rule check_md5_short_read:
    input:
        md5_check_json = os.path.join('../00.raw_data',config['short_read_convert_md5'],"raw_data_md5.json"),
        link_r1 = expand(os.path.join('../00.raw_data',
                                      config['short_read_convert_md5'],
                                      "{sample}/{sample}_R1.fq.gz"),
                                      sample=samples.keys()),
        link_r2 = expand(os.path.join('../00.raw_data',
                                      config['short_read_convert_md5'],
                                      "{sample}/{sample}_R2.fq.gz"),
                                      sample=samples.keys()),
    output:
        md5_check = "../01.qc/md5_check_short_read.tsv",
    message:
        "Running md5 check on raw data files on {input.md5_check_json}",
    benchmark:
        "../benchmarks/md5_check_short_read_benchmark.txt",
    log:
        "../logs/01.qc/md5_check_short_read.log",
    params:
        md5_check = os.path.join('../00.raw_data',config['short_read_convert_md5']),
        log_file = "../logs/01.qc/md5_check_short_read.log",
    threads: 
        config['threads']['short_md5_check'],
    shell:
        """
        ./scripts/json_md5_verifier/target/release/json_md5_verifier -t  {threads} \
                -i {input.md5_check_json} \
                -b {params.md5_check} \
                -o {output.md5_check} \
                --log-file {params.log_file} &> {log}
        """
# ------------- long read check ------------ #
if config['long_read']:
    rule seq_preprocessor_short_read:
        input:
            md5 = config['long_read_raw_data_path'],
            # md5 = get_all_input_dirs(load_samples.keys(),config = config,raw_data_path_config = 'long_read_raw_data_path'),
        output:
            md5_check = directory(os.path.join('../00.raw_data',config['long_read_convert_md5'])),
            md5_check_json = os.path.join('../00.raw_data',config['long_read_convert_md5'],
                                           "long_raw_data_md5.json"),
            link_se = expand(os.path.join('../00.raw_data',
                                        config['short_read_convert_md5'],
                                        "{long_sample}/{long_sample}.fq.gz"),
                                        long_sample = long_read_samples.keys()),
        message:
            "Running seq_preprocessor on raw data data",
        benchmark:
            "../benchmarks/seq_preprocessor.txt",
        params:
            raw_data_path = config['long_read_raw_data_path'],
            md5 = config['long_read_md5'],
        log:
            "../logs/01.qc/seq_preprocessor.txt",
        threads: 1
        shell:
            """
            ./scripts/seq_preprocessor/target/release/seq_preprocessor -i  {params.raw_data_path} \
                    -o {output.md5_check} \
                    --library-type long-read \
                    --md5-name {params.md5} \
                    --json-report {output.md5_check_json} &> {log}
            """
    rule check_md5_short_read:
        input:
            md5_check_json = os.path.join('../00.raw_data',config['short_read_convert_md5'],"long_raw_data_md5.json"),
            link_se = expand(os.path.join('../00.raw_data',
                                        config['short_read_convert_md5'],
                                        "{long_sample}/{long_sample}.fq.gz"),
                                        long_sample=long_read_samples.keys()),
        output:
            md5_check = "../01.qc/md5_check_short_read.tsv",
        message:
            "Running md5 check on raw data files on {input.md5_check_json}",
        benchmark:
            "../benchmarks/md5_check_short_read_benchmark.txt",
        log:
            "../logs/01.qc/md5_check_short_read.log",
        params:
            md5_check = os.path.join('../00.raw_data',config['short_read_convert_md5']),
            log_file = "../logs/01.qc/md5_check_short_read.log",
        threads: 
            config['threads']['md5_check'],
        shell:
            """
            ./scripts/json_md5_verifier/target/release/json_md5_verifier -t  {threads} \
                    -i {input.md5_check_json} \
                    -b {params.md5_check} \
                    -o {output.md5_check} \
                    --log-file {params.log_file} &> {log}
            """

# ------- rule ------- #