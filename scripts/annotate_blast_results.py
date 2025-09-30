#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from loguru import logger
import pandas as pd
import argparse
import sys


# --- Loguru 配置 ---
# 移除默认的 handler，添加一个更美观、信息更丰富的 handler
logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>"
)
# --------------------

def get_best_hits(blast_result_path):
    """
    从 BLAST/DIAMOND 的比对结果中提取每个 Query 的最佳匹配项。
    假定输入文件已经按照 E-value 和 bit score 排序，因此每个 Query 的第一个出现就是最佳匹配。
    """
    best_hits = {}
    
    logger.info(f"开始解析比对文件: {blast_result_path}")
    
    with open(blast_result_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            
            query_id = fields[0]
            
            if query_id not in best_hits:
                subject_id = fields[1]
                clean_subject_id = subject_id.split('.')[0]
                
                hit_data = {
                    'query_id': query_id,
                    'subject_id': subject_id,
                    'clean_subject_id': clean_subject_id,
                    'identity': float(fields[2]),
                    'e_value': float(fields[10]),
                    'bit_score': float(fields[11]),
                }
                best_hits[query_id] = hit_data

    logger.success(f"解析完成, 找到 {len(best_hits)} 个独立查询序列的最佳匹配。")
    return pd.DataFrame(list(best_hits.values()))

def main():
    parser = argparse.ArgumentParser(description="根据 DIAMOND/BLAST 比对结果和 UniProt 注释文件，生成最终的功能注释表。")
    parser.add_argument('--blast_result', required=True, help="DIAMOND/BLAST 的 tabular 输出文件 (outfmt 6)。")
    parser.add_argument('--uniprot_map', required=True, help="从 UniProt 下载的 TSV 格式注释文件。")
    parser.add_argument('--output', required=True, help="最终输出的注释文件名 (TSV 格式)。")
    
    args = parser.parse_args()

    try:
        # 1. 读取 UniProt 注释文件
        logger.info(f"正在加载 UniProt 注释文件: {args.uniprot_map}")

        # --- DtypeWarning 解决方案 ---
        # 明确指定可能引起警告的列的数据类型为字符串 (str)
        dtype_spec = {
            'Biotechnological use': str,
            'Pharmaceutical use': str,
            'Mutagenesis': str,
            'Involvement in disease': str,
            'Toxic dose': str
        }
        uniprot_df = pd.read_csv(args.uniprot_map, sep='\t', dtype=dtype_spec)
        # ----------------------------

        uniprot_df.rename(columns={'Entry': 'clean_subject_id'}, inplace=True)
        uniprot_df.set_index('clean_subject_id', inplace=True)
        logger.success("UniProt 注释文件加载成功。")

        # 2. 从比对结果中提取最佳匹配
        best_hits_df = get_best_hits(args.blast_result)
        
        # 3. 将最佳匹配结果与 UniProt 注释信息合并
        logger.info("正在合并比对结果和注释信息...")
        final_df = best_hits_df.merge(uniprot_df, on='clean_subject_id', how='left')
        
        # 4. 保存最终结果
        logger.info(f"正在保存最终注释结果到: {args.output}")
        
        output_columns = [
            'query_id', 'subject_id', 'identity', 'e_value', 'bit_score', 'Reviewed',
            'Entry Name', 'Protein names', 'Gene Names', 'Organism', 'Length',
            'Gene Ontology (biological process)', 'Gene Ontology (cellular component)',
            'Gene Ontology (molecular function)', 'Allergenic Properties', 
            'Biotechnological use', 'Disruption phenotype', 'Pharmaceutical use',
            'Mutagenesis', 'Involvement in disease', 'Toxic dose'
        ]
        
        existing_columns = [col for col in output_columns if col in final_df.columns]
        
        final_df[existing_columns].to_csv(args.output, sep='\t', index=False, na_rep='NA')
        logger.success("全部完成！")

    except FileNotFoundError as e:
        logger.error(f"文件未找到: {e.filename}")
        sys.exit(1)
    except KeyError as e:
        logger.error(f"注释文件中未找到预期的列: {e}。请检查文件是否正确。")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"脚本运行中发生未知错误: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()