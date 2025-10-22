#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os

def clean_kaas_annotations(input_file, output_file):
    """
    读取KAAS的原始注释文件，
    过滤掉没有KO注释的行，
    生成一个干净的 'gene2ko' 映射文件。
    """
    print(f"--- 开始处理KAAS注释文件: {input_file} ---")
    
    processed_count = 0
    skipped_count = 0
    
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            
            # 为输出文件写入表头，这在下游分析中很有用
            outfile.write("TranscriptID\tKO_ID\n")
            
            for line in infile:
                # 去除行首尾的空白（如换行符）
                line_stripped = line.strip()
                
                # 跳过空行
                if not line_stripped:
                    continue
                
                # 尝试按制表符（\t）分割
                parts = line_stripped.split('\t')
                
                if len(parts) == 2:
                    # 长度为2，表示是 "ID \t KO" 的标准格式
                    gene_id = parts[0]
                    ko_id = parts[1].strip() # 再次strip以防万一
                    
                    # 关键：检查KO ID是否存在且不为空
                    if ko_id:
                        outfile.write(f"{gene_id}\t{ko_id}\n")
                        processed_count += 1
                    else:
                        # KO ID 为空 (e.g., "rb_1.p1\t")
                        skipped_count += 1
                        
                elif len(parts) == 1:
                    # 长度为1，表示只有ID，没有KO (e.g., "rb_1.p1")
                    skipped_count += 1
                
                else:
                    # 格式奇怪的行
                    print(f"警告: 跳过格式不规范的行: {line_stripped}", file=sys.stderr)

    except FileNotFoundError:
        print(f"错误: 找不到输入文件 '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"处理过程中发生错误: {e}", file=sys.stderr)
        sys.exit(1)

    print("\n--- 处理完成 ---")
    print(f"成功写入的注释条目: {processed_count}")
    print(f"跳过的未注释条目: {skipped_count}")
    print(f"下游分析文件已保存至: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="清洗KAAS的KEGG注释结果，为下游分析做准备。")
    
    parser.add_argument("-i", "--input", 
                        required=True, 
                        help="KAAS输出的原始注释文件 (例如: kaas_results.txt)")
    
    parser.add_argument("-o", "--output", 
                        default="gene2ko.tsv", 
                        help="清洗后用于下游分析的 'gene-to-KO' 文件名 (默认: gene2ko.tsv)")
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input):
        print(f"错误: 输入文件 {args.input} 不存在。", file=sys.stderr)
        sys.exit(1)
        
    clean_kaas_annotations(args.input, args.output)