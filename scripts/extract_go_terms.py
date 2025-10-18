import pandas as pd
import argparse
import logging
import re
import sys

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
    # 例如: "diterpenoid biosynthetic process [GO:0016102]"
    # (.*) 会捕获描述部分, (GO:\d+) 会捕获GO ID部分
    go_pattern = re.compile(r'(.*) \[(GO:\d+)\]')
    
    try:
        logger.info(f"正在读取输入文件: {args.input}")
        # 使用pandas读取制表符分隔的文件
        df = pd.read_csv(args.input, sep='\t')
        logger.info(f"成功加载 {len(df)} 条记录。")
    except FileNotFoundError:
        logger.error(f"错误：输入文件未找到 -> {args.input}")
        sys.exit(1) # 退出脚本
    except Exception as e:
        logger.error(f"读取文件时发生错误: {e}")
        sys.exit(1)

    # 创建一个空列表来存储所有提取出的GO信息
    extracted_data = []

    logger.info("开始提取 GO 条目...")
    # 遍历DataFrame的每一行
    for index, row in df.iterrows():
        gene_id = row['query_id']
        
        # 遍历三个GO列
        for col_name, go_type in go_columns.items():
            # 检查单元格是否为空 (pd.notna 处理 NaN 等空值)
            if col_name in row and pd.notna(row[col_name]):
                # 获取单元格内容并按分号切分，得到多个GO条目
                entries = str(row[col_name]).split(';')
                
                for entry in entries:
                    entry = entry.strip() # 去除首尾空格
                    if not entry:
                        continue # 跳过空条目
                    
                    # 使用正则表达式进行匹配
                    match = go_pattern.match(entry)
                    if match:
                        # 如果匹配成功，提取描述和GO ID
                        description, go_id = match.groups()
                        # 将结果添加到列表中
                        extracted_data.append([
                            gene_id,
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
        # 将结果保存到指定的输出文件，格式为TSV，不包含索引
        logger.info(f"正在将结果保存到: {args.output}")
        output_df.to_csv(args.output, sep='\t', index=False)
        logger.info("处理完成！")
    except Exception as e:
        logger.error(f"保存文件时发生错误: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()