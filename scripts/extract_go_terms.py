#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import argparse
import logging
import pandas as pd

def setup_logger():
    """配置日志记录器，用于在控制台输出脚本运行状态。"""
    # 设置日志记录的基本配置
    logging.basicConfig(
        level=logging.INFO,  # 日志级别为INFO
        stream=sys.stdout,   # 输出到标准输出（控制台）
        format='%(asctime)s - %(levelname)s - %(message)s', # 日志格式
        datefmt='%Y-%m-%d %H:%M:%S' # 时间格式
    )
    return logging.getLogger(__name__)

def parse_arguments():
    """解析命令行参数。"""
    parser = argparse.ArgumentParser(
        description="从 UniProt 注释结果中提取 GO (Gene Ontology) 富集信息。"
    )
    # 定义输入文件参数
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="输入的 UniProt 注释文件路径 (TSV 格式)。"
    )
    # 定义输出文件参数
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="输出提取后的 GO 信息文件路径 (TSV 格式)。"
    )
    return parser.parse_args()

def main():
    """主执行函数。"""
    # 初始化
    args = parse_arguments()
    logger = setup_logger()
    
    logger.info("脚本开始运行...")
    
    # 定义GO注释所在的列名和它们对应的GO类型
    go_columns = {
        'Gene Ontology (biological process)': 'biological_process',
        'Gene Ontology (cellular component)': 'cellular_component',
        'Gene Ontology (molecular function)': 'molecular_function'
    }
    
    # 用于从字符串中匹配 "描述 [GO:ID]" 格式的正则表达式
    go_pattern = re.compile(r'(.*) \[(GO:\d+)\]')
    
    try:
        logger.info(f"正在读取输入文件: {args.input}")
        df = pd.read_csv(args.input, sep='\t',low_memory = False)
        logger.info(f"成功加载 {len(df)} 条记录。")
    except FileNotFoundError:
        logger.error(f"错误：输入文件未找到 -> {args.input}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"读取文件时发生错误: {e}")
        sys.exit(1)

    # 检查 'Entry Name' 列是否存在
    if 'Entry Name' not in df.columns:
        logger.error(f"错误: 输入文件中未找到 'Entry Name' 列。")
        sys.exit(1)
    if 'query_id' not in df.columns:
        logger.error(f"错误: 输入文件中未找到 'query_id' 列。")
        sys.exit(1)

    extracted_data = []

    logger.info("开始提取 GO 条目...")
    # 遍历DataFrame的每一行
    for index, row in df.iterrows():
        gene_id = row['query_id']
        
        # --- 修正点 (开始) ---
        # 1. 从 row 中获取 'Entry Name' 的值
        entry_name_val = row['Entry Name']
        
        # 2. 检查 'Entry Name' 是否为空 (NaN)
        if pd.notna(entry_name_val):
            # 如果不为空，拼接成 "gene_id"_"Entry_Name" 的形式
            combined_id = f"{str(entry_name_val)}_{gene_id}"
        else:
            # 如果 'Entry Name' 为空，则只使用 gene_id 作为ID
            combined_id = gene_id
        # --- 修正点 (结束) ---
            
        # 遍历三个GO列
        for col_name, go_type in go_columns.items():
            if col_name in row and pd.notna(row[col_name]):
                entries = str(row[col_name]).split(';')
                
                for entry in entries:
                    entry = entry.strip()
                    if not entry:
                        continue
                    
                    match = go_pattern.match(entry)
                    if match:
                        description, go_id = match.groups()
                        
                        # --- 修正点：使用我们新创建的 combined_id ---
                        extracted_data.append([
                            combined_id,  # <-- 使用修正后的ID
                            go_type,
                            go_id,
                            description.strip()
                        ])

    if not extracted_data:
        logger.warning("未找到任何 GO 条目。输出文件将为空。")
    else:
        logger.info(f"成功提取 {len(extracted_data)} 条 GO 记录。")

    # 将结果列表转换为一个新的DataFrame
    output_df = pd.DataFrame(
        extracted_data,
        columns=['GeneID', 'GO_Type', 'GO_ID', 'GO_Description']
    )
    
    try:
        logger.info(f"正在将结果保存到: {args.output}")
        output_df.to_csv(args.output, sep='\t', index=False)
        logger.info("处理完成！")
    except Exception as e:
        logger.error(f"保存文件时发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()