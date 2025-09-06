#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import pandas as pd

if len(sys.argv) != 4:
    print("User: python3 E90_filter.py <quant.sf_path> <percentage> <output_ids_file>")
    print("Example: python filter_e90.py salmon_output/quant.sf 0.90 core_transcript_ids.txt")
    sys.exit(1)
# --- 参数读取 ---
quant_file = sys.argv[1]
percent_threshold = float(sys.argv[2]) 
output_file = sys.argv[3]

# --- 核心逻辑 ---
print(f"正在读取定量文件: {quant_file}...")
# 读取 Salmon 的 quant.sf 文件
df = pd.read_csv(quant_file, sep='\t')

# 按 TPM 值从高到低排序
df_sorted = df.sort_values(by='TPM', ascending=False)

# 计算总 TPM
total_tpm = df_sorted['TPM'].sum()
print(f"总 TPM 值为: {total_tpm:.2f}")

# 计算阈值
tpm_cutoff = total_tpm * percent_threshold
print(f"将选取转录本直到 TPM 累积和达到: {tpm_cutoff:.2f} (总和的 {percent_threshold * 100}%)")

# 循环以确定核心转录本
cumulative_tpm = 0
core_transcripts_ids = []

for index, row in df_sorted.iterrows():
    if cumulative_tpm < tpm_cutoff:
        cumulative_tpm += row['TPM']
        core_transcripts_ids.append(row['Name'])
    else:
        # 一旦达到阈值，就停止添加
        break

print(f"已选取 {len(core_transcripts_ids)} 个转录本，其累积 TPM 为 {cumulative_tpm:.2f}")

# 将选中的转录本 ID 写入文件
with open(output_file, 'w') as f:
    for transcript_id in core_transcripts_ids:
        f.write(f"{transcript_id}\n")

print(f"核心转录本 ID 已保存至: {output_file}")