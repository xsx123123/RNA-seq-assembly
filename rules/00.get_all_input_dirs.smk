#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# loading packages
from loguru import logger
from pathlib import Path
from typing import Dict, Union
from rich import print as rich_print
# ----------------------------- #
def get_sample_data_dir(sample_id:str = None,
                        config:dict = None,
                        raw_data_path_config:str = "raw_data_path") -> str:
    """
    根据 *具体的* sample_id (e.g., "Sample_A"),
    查找 *包含* fastq 文件的 *目录*。
    
    注意：我修改了它，使其不再依赖 wildcards，
    而是直接接收 sample_id 字符串。
    """
    
    for base_dir in config[raw_data_path_config]:
        sample_dir = os.path.join(base_dir, sample_id)
        if os.path.isdir(sample_dir):
            return sample_dir
                
    raise FileNotFoundError(f"无法在 {config[raw_data_path_config]} 中找到 {sample_id} ({sample_path_component}) 的数据目录")

def get_all_input_dirs(sample_keys:str = None,
                       config:dict = config) -> list:
    """
    遍历所有样本 ID，调用 get_sample_data_dir，
    返回一个包含所有数据目录的列表。
    """
    dir_list = []
    for sample_id in sample_keys:
        dir_list.append(get_sample_data_dir(sample_id,config = config))
    print(dir_list)
    return list(set(dir_list))