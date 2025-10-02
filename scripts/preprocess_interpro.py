#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys
import argparse
from loguru import logger
from tqdm import tqdm
from pandarallel import pandarallel

def main():
    """
    主函数，用于解析命令行参数并执行处理流程。
    """
    # --- 配置 loguru 日志 ---
    logger.remove()
    logger.add(
        sys.stderr, 
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}"
    )

    # --- 设置命令行参数解析 ---
    parser = argparse.ArgumentParser(
        description="Preprocess InterProScan TSV output to generate a summary table and a GO mapping file for enrichment analysis."
    )
    parser.add_argument("input_file", help="Path to the InterProScan TSV file.")
    parser.add_argument(
        "--summary_out", 
        default="protein_annotation_summary.tsv", 
        help="Output file name for the summary table (default: protein_annotation_summary.tsv)."
    )
    parser.add_argument(
        "--go_map_out", 
        default="protein_to_go_mapping.tsv", 
        help="Output file name for the GO enrichment mapping file (default: protein_to_go_mapping.tsv)."
    )
    # 添加一个用于控制CPU核心数的参数
    parser.add_argument(
        "--cores", 
        type=int, 
        default=4, 
        help="Number of CPU cores to use for parallel processing (default: 4)."
    )
    args = parser.parse_args()

    # --- 2. 初始化 pandarallel ---
    # nb_workers 参数指定了使用的 CPU 核心数，并显示进度条
    logger.info(f"Initializing parallel processing with {args.cores} CPU cores.")
    pandarallel.initialize(nb_workers=args.cores, progress_bar=True)

    # --- 读取并预处理数据 ---
    logger.info(f"Reading InterProScan file: {args.input_file}")
    logger.info("This may take a while for large files...")
    
    column_names = [
        'Protein_ID', 'MD5', 'Length', 'Database', 'DB_Accession',
        'Description', 'Start', 'Stop', 'E-Value', 'Status', 'Date',
        'InterPro_Accession', 'InterPro_Description', 'GO_Terms', 'Pathways'
    ]
    
    try:
        df = pd.read_csv(
            args.input_file, 
            sep='\t', 
            header=None, 
            names=column_names,
            on_bad_lines='warn',
            low_memory=False
        )
        logger.success("File read successfully.")
    except FileNotFoundError:
        logger.error(f"Input file not found: {args.input_file}")
        sys.exit(1)

    # --- 生成富集分析使用的表 (这个操作通常很快，无需并行) ---
    logger.info("Generating GO mapping file for enrichment analysis...")
    
    go_df = df[df['GO_Terms'] != '-'].copy()
    go_df['GO_Terms'] = go_df['GO_Terms'].str.split('|')
    enrichment_table = go_df.explode('GO_Terms')
    enrichment_table = enrichment_table[['Protein_ID', 'GO_Terms']].drop_duplicates().reset_index(drop=True)
    
    enrichment_table.to_csv(args.go_map_out, sep='\t', index=False)
    logger.success(f"GO mapping file saved to: {args.go_map_out}")

    # --- 生成总注释表 (并行处理) ---
    logger.info("Generating summary annotation table using parallel processing...")

    def summarize_protein(group):
        def get_terms(db_name):
            subset = group[group['Database'] == db_name]
            if not subset.empty:
                return ';'.join(
                    f"{acc}({desc})" for acc, desc in subset[['DB_Accession', 'Description']].drop_duplicates().itertuples(index=False)
                )
            return '-'

        all_gos = set()
        group['GO_Terms'].dropna().str.split('|').apply(lambda x: all_gos.update(x))
        all_gos.discard('-')

        interpro_terms = ';'.join(group['InterPro_Accession'].dropna().unique()).replace('-', '')

        return pd.Series({
            'Length': group['Length'].iloc[0],
            'Pfam_Annotation': get_terms('Pfam'),
            'PANTHER_Annotation': get_terms('PANTHER'),
            'InterPro_Entries': interpro_terms if interpro_terms else '-',
            'GO_Terms': ';'.join(sorted(list(all_gos))) if all_gos else '-'
        })

    # --- 3. 使用 .parallel_apply() 代替 .progress_apply() ---
    summary_table = df.groupby('Protein_ID').parallel_apply(summarize_protein)
    
    summary_table.to_csv(args.summary_out, sep='\t')
    logger.success(f"Summary table saved to: {args.summary_out}")
    logger.info("All tasks completed.")


if __name__ == "__main__":
    main()
