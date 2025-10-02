import pandas as pd

# --- 1. 文件路径配置 ---
# 请根据你的实际文件名修改这里的路径
file_paths = {
    # 包含所有转录本ID的文件 (例如，从转录组fasta文件提取的ID列表)
    "transcript_ids": "transcript_ids.txt",
    
    # BLASTp/BLASTX 的输出结果 (推荐使用 outfmt 6 格式)
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    "blast": "blastp_results.tsv",
    
    # hmmscan 的输出结果 (必须是 --domtblout 格式)
    "pfam": "pfam_results.domtblout",
    
    # GO 注释文件 (假设是两列：Transcript_ID GO_ID)
    "go": "go_annotations.tsv"
}

output_file = "final_annotations.tsv"

# --- 2. 解析函数定义 ---

def load_all_transcripts(filepath):
    """从一个单列表文件中加载所有转录本ID，作为基础DataFrame。"""
    print(f"Loading all transcript IDs from {filepath}...")
    # 假设文件每行一个ID，没有表头
    df = pd.read_csv(filepath, header=None, names=['transcript_id'])
    print(f"Found {len(df)} total transcripts.")
    return df

def parse_blast(filepath):
    """解析BLAST outfmt 6结果，并为每个query保留最佳hit。"""
    print(f"Parsing BLAST results from {filepath}...")
    try:
        # 定义BLAST outfmt 6的标准列名
        col_names = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
        ]
        df = pd.read_csv(filepath, sep='\t', header=None, names=col_names)
        
        # 按bitscore降序、evalue升序排序
        df_sorted = df.sort_values(by=['bitscore', 'evalue'], ascending=[False, True])
        
        # 去除重复的qseqid，只保留第一个（即最佳hit）
        best_hits = df_sorted.drop_duplicates('qseqid', keep='first')
        
        # 选择我们感兴趣的列
        best_hits = best_hits[['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']]
        best_hits = best_hits.rename(columns={'qseqid': 'transcript_id', 'sseqid': 'blast_hit'})
        print(f"Processed {len(best_hits)} unique BLAST hits.")
        return best_hits
    except FileNotFoundError:
        print(f"Warning: BLAST file not found at {filepath}. Skipping.")
        return pd.DataFrame({'transcript_id': []})


def parse_pfam(filepath):
    """解析hmmscan --domtblout 结果，合并一个转录本的多个domain。"""
    print(f"Parsing Pfam results from {filepath}...")
    try:
        # hmmscan的domtblout输出是空格分隔的，且前面有注释行(#)
        with open(filepath, 'r') as f:
            lines = [line for line in f if not line.startswith('#')]
        
        if not lines:
            print("Warning: Pfam file is empty or contains no data rows. Skipping.")
            return pd.DataFrame({'transcript_id': []})

        # 将清理过的行读入DataFrame
        col_names = [
            'domain_name', 'domain_accession', 'tlen', 'query_name', 'query_accession', 'qlen',
            'evalue', 'score', 'bias', '#', 'of', 'c-evalue', 'i-evalue', 'score_dom', 'bias_dom',
            'from_hmm', 'to_hmm', 'from_ali', 'to_ali', 'from_env', 'to_env', 'acc', 'description'
        ]
        # 使用正则表达式匹配一个或多个空格作为分隔符
        df = pd.read_csv(pd.io.common.StringIO(''.join(lines)), sep='\s+', header=None, names=col_names)
        
        # 筛选可靠的hit (e.g., E-value < 1e-5)
        df_filtered = df[df['evalue'].astype(float) < 1e-5]

        # 按query_name分组，并将domain_name合并成一个字符串
        # df_agg = df_filtered.groupby('query_name')['domain_name'].apply(lambda x: ';'.join(x)).reset_index()
        # 同时保留domain name和对应的e-value
        df_filtered['domain_info'] = df_filtered['domain_name'] + '(' + df_filtered['evalue'].astype(str) + ')'
        df_agg = df_filtered.groupby('query_name')['domain_info'].apply(';'.join).reset_index()
        
        df_agg = df_agg.rename(columns={'query_name': 'transcript_id', 'domain_info': 'pfam_domains'})
        print(f"Processed Pfam annotations for {len(df_agg)} transcripts.")
        return df_agg
    except FileNotFoundError:
        print(f"Warning: Pfam file not found at {filepath}. Skipping.")
        return pd.DataFrame({'transcript_id': []})

def parse_go(filepath):
    """解析GO注释文件，合并一个转录本的多个GO term。"""
    print(f"Parsing GO results from {filepath}...")
    try:
        # 假设文件是两列: transcript_id, go_term，以tab分隔
        df = pd.read_csv(filepath, sep='\t', header=None, names=['transcript_id', 'go_term'])
        
        # 按transcript_id分组，并将go_term合并
        df_agg = df.groupby('transcript_id')['go_term'].apply(';'.join).reset_index()
        print(f"Processed GO annotations for {len(df_agg)} transcripts.")
        return df_agg
    except FileNotFoundError:
        print(f"Warning: GO file not found at {filepath}. Skipping.")
        return pd.DataFrame({'transcript_id': []})

# --- 3. 主逻辑 ---

def main():
    """主函数，执行所有解析和合并步骤。"""
    # 加载所有转录本ID作为基础
    base_df = load_all_transcripts(file_paths['transcript_ids'])
    
    # 解析各个注释文件
    blast_df = parse_blast(file_paths['blast'])
    pfam_df = parse_pfam(file_paths['pfam'])
    go_df = parse_go(file_paths['go'])
    
    # 开始合并
    print("\nMerging all annotation data...")
    # 使用左合并(left merge)，确保所有原始转录本都保留
    merged_df = pd.merge(base_df, blast_df, on='transcript_id', how='left')
    merged_df = pd.merge(merged_df, pfam_df, on='transcript_id', how='left')
    merged_df = pd.merge(merged_df, go_df, on='transcript_id', how='left')
    
    # 将所有NaN (Not a Number) 的缺失值替换为 '.'
    merged_df = merged_df.fillna('.')
    
    # 保存结果
    print(f"Writing final merged annotations to {output_file}...")
    merged_df.to_csv(output_file, sep='\t', index=False)
    
    print("\nDone! ✨")

if __name__ == "__main__":
    main()