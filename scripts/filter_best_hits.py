#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import sys
import argparse
import logging

def setup_logging():
    """配置日志记录"""
    logging.basicConfig(
        level=logging.INFO,  # 设置日志级别为 INFO
        format='%(asctime)s [%(levelname)s] - %(message)s', # 日志格式
        datefmt='%Y-%m-%d %H:%M:%S' # 日期时间格式
    )

def find_best_hits(input_file, output_file):
    """
    读取 BLAST/UniProt 的注释结果，
    并根据 'bit_score' 为每个转录本筛选出最佳的匹配。
    """
    logging.info(f"正在读取文件: {input_file} ...")
    try:
        # 读取TSV文件
        df = pd.read_csv(input_file, sep='\t')
    except FileNotFoundError:
        logging.error(f"错误：找不到文件 {input_file}")
        logging.error("请确保提供了正确的文件路径。")
        sys.exit(1)
    except Exception as e:
        logging.error(f"读取文件时发生未知错误: {e}")
        sys.exit(1)

    logging.info(f"总共加载了 {len(df)} 条注释记录。")

    # 确保 'bit_score' 列是数值类型
    if 'bit_score' not in df.columns:
        logging.error("错误：输入文件中未找到 'bit_score' 列。")
        logging.error("请确保文件包含 bit_score 用于排序。")
        sys.exit(1)
        
    df['bit_score'] = pd.to_numeric(df['bit_score'], errors='coerce')
    df = df.dropna(subset=['bit_score']) # 移除 bit_score 为空的行

    logging.info("正在从 'query_id' 提取转录本ID (例如 'rb_1.p1' -> 'rb_1')...")
    
    # 从 'rb_1.p1' 中提取 'rb_1'
    df['transcript_id'] = df['query_id'].str.rsplit('.', n=1).str[0]

    logging.info("正在按转录本ID分组，并查找 'bit_score' 最高的匹配... 🥇")
    
    # 核心逻辑：
    # 1. 按 'transcript_id' 分组
    # 2. 在每个组内找到 'bit_score' 最大值的索引
    try:
        best_indices = df.groupby('transcript_id')['bit_score'].idxmax()
    except Exception as e:
        logging.error(f"根据 'bit_score' 分组排序时出错: {e}")
        logging.error("请检查 'transcript_id' 和 'bit_score' 列的数据。")
        sys.exit(1)

    # 3. 使用索引选出“最佳”行
    best_hits_df = df.loc[best_indices]

    # 将 'transcript_id' 列移动到最前面，方便查看
    cols = list(best_hits_df.columns)
    cols.insert(0, cols.pop(cols.index('transcript_id')))
    best_hits_df = best_hits_df[cols]

    logging.info(f"筛选完毕！共找到 {len(best_hits_df)} 个转录本的最佳注释。")

    # 保存结果
    try:
        best_hits_df.to_csv(output_file, sep='\t', index=False, na_rep='NA')
        logging.info(f"🎉 成功！最佳注释已保存到: {output_file}")
    except Exception as e:
        logging.error(f"保存文件到 {output_file} 时出错: {e}")
        sys.exit(1)

# --- 脚本主入口 ---
if __name__ == "__main__":
    # 1. 配置日志
    setup_logging()

    # 2. 配置参数解析
    parser = argparse.ArgumentParser(
        description="根据 'bit_score' 为每个转录本筛选最佳的 UniProt 注释。",
        epilog="示例: python filter_best_hits.py -i TD2_pep_matches_annotated.tsv -o best_hits.tsv"
    )
    
    # 添加 -i (输入) 参数
    parser.add_argument(
        '-i', '--input',
        required=True,  # 设为必需参数
        help="输入的TSV注释文件 (例如 Trinotate 的 BLASTx/UniProt 报告)"
    )
    
    # 添加 -o (输出) 参数
    parser.add_argument(
        '-o', '--output',
        required=True,  # 设为必需参数
        help="筛选后输出的最佳注释TSV文件"
    )

    # 解析命令行参数
    args = parser.parse_args()

    # 3. 运行主函数
    find_best_hits(args.input, args.output)