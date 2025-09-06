#!/usr/bin/env python
import pandas as pd
import argparse
from pathlib import Path
import sys

def merge_salmon_quants(file_paths, output_file):
    """
    Merges TPM values from multiple Salmon quant.sf files into a single matrix.

    Args:
        file_paths (list): A list of paths to quant.sf files.
        output_file (str): The path for the output matrix file.
    """
    all_dataframes = []

    if not file_paths:
        print("错误：没有提供任何输入文件。", file=sys.stderr)
        sys.exit(1)

    print(f"开始处理 {len(file_paths)} 个样本文件...")

    for f_path in file_paths:
        try:
            # 使用 pathlib 来方便地获取父目录名作为样本名
            # 这与 Trinity 脚本的 --name_sample_by_basedir 行为一致
            p = Path(f_path)
            sample_name = p.parent.name
            
            # 读取 quant.sf 文件，它是一个制表符分隔的文件
            df = pd.read_csv(f_path, sep='\t', usecols=['Name', 'TPM'])
            
            # 将 'TPM' 列重命名为样本名，以便合并
            df.rename(columns={'TPM': sample_name}, inplace=True)
            
            # 将 'Name' (转录本ID) 列设为索引，这是合并的关键
            df.set_index('Name', inplace=True)
            
            all_dataframes.append(df)
            print(f"  - 已处理样本: {sample_name}")

        except FileNotFoundError:
            print(f"警告：找不到文件 {f_path}，已跳过。", file=sys.stderr)
        except Exception as e:
            print(f"处理文件 {f_path} 时出错: {e}", file=sys.stderr)

    if not all_dataframes:
        print("错误：没有成功处理任何文件，无法生成矩阵。", file=sys.stderr)
        sys.exit(1)

    # 使用 concat 函数沿着列（axis=1）合并所有的 DataFrame
    # Pandas 会自动根据索引（转录本ID）对齐数据
    print("\n正在合并所有样本的数据...")
    merged_df = pd.concat(all_dataframes, axis=1)

    # 将合并过程中产生的 NaN 值（代表某个样本中没有该转录本的表达）填充为 0
    merged_df.fillna(0, inplace=True)

    # 将合并后的矩阵写入输出文件，以制表符分隔
    print(f"正在将结果写入: {output_file}")
    merged_df.to_csv(output_file, sep='\t')
    
    print("处理完成！")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="将多个 Salmon quant.sf 文件的 TPM 值合并成一个表达矩阵。"
    )
    parser.add_argument(
        'infiles', 
        nargs='+',  # '+' 表示接受一个或多个位置参数
        help="一个或多个 quant.sf 文件的路径。"
    )
    parser.add_argument(
        '-o', '--output', 
        required=True, 
        help="输出的矩阵文件名。"
    )
    
    args = parser.parse_args()
    
    merge_salmon_quants(args.infiles, args.output)