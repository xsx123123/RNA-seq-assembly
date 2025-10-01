#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import pandas as pd

# 检查命令行参数数量是否正确
if len(sys.argv) != 4:
    print("用法: python3 E90_filter.py <tpm_matrix_path> <percentage> <output_ids_file>")
    print("示例: python3 E90_filter.py all_samples.TPM.matrix 0.90 core_transcript_ids.txt")
    sys.exit(1)

# --- 参数读取 ---
matrix_file = sys.argv[1]
percent_threshold = float(sys.argv[2]) 
output_file = sys.argv[3]

# --- 核心逻辑 ---
print(f"正在读取表达矩阵文件: {matrix_file}...")
try:
    # 读取 TPM 矩阵文件
    df = pd.read_csv(matrix_file, sep='\t')
except FileNotFoundError:
    print(f"错误: 文件 '{matrix_file}' 未找到。", file=sys.stderr)
    sys.exit(1)

# 将第一列（转录本ID）设为索引
df.set_index(df.columns[0], inplace=True)

# 计算每个转录本在所有样本中的平均 TPM
print("正在计算每个转录本的平均 TPM...")
df['Mean_TPM'] = df.mean(axis=1)

# 创建一个新的 DataFrame 用于排序，只保留 Name 和 Mean_TPM
mean_tpm_df = df[['Mean_TPM']].reset_index()
mean_tpm_df.rename(columns={mean_tpm_df.columns[0]: 'Name'}, inplace=True)


# 按 Mean_TPM 值从高到低排序
df_sorted = mean_tpm_df.sort_values(by='Mean_TPM', ascending=False)

# 计算总的平均 TPM 之和
total_mean_tpm = df_sorted['Mean_TPM'].sum()
print(f"所有转录本的平均 TPM 总值为: {total_mean_tpm:.2f}")

# 计算筛选阈值
tpm_cutoff = total_mean_tpm * percent_threshold
print(f"将选取转录本直到平均 TPM 累积和达到: {tpm_cutoff:.2f} (总和的 {percent_threshold * 100:.1f}%)")

# 循环以确定核心转录本
cumulative_tpm = 0
core_transcripts_ids = []

for index, row in df_sorted.iterrows():
    # 如果累积值小于阈值，则继续添加
    if cumulative_tpm < tpm_cutoff:
        cumulative_tpm += row['Mean_TPM']
        core_transcripts_ids.append(row['Name'])
    else:
        # 一旦达到或超过阈值，就停止添加
        break

print(f"已选取 {len(core_transcripts_ids)} 个转录本，其累积平均 TPM 为 {cumulative_tpm:.2f}")

# 将选中的转录本 ID 写入文件
with open(output_file, 'w') as f:
    for transcript_id in core_transcripts_ids:
        f.write(f"{transcript_id}\n")

print(f"核心转录本 ID 已保存至: {output_file}")