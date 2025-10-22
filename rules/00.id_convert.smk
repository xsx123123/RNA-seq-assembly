#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from loguru import logger
from rich import inspect,print
# Load the samples CSV file into a DataFrame
samples_df = pd.read_csv(config["sample_csv"])
logger.info(f"Loaded samples from {config["sample_csv"]}")
# Convert the DataFrame to a dictionary for easy access
load_samples = samples_df.set_index("sample").to_dict(orient="index")
logger.info(f"Analyzed {len(load_samples)} short-read samples from the CSV file.")
# ------------- #
long_read_samples_df = pd.read_csv(config["long_read_sample_csv"])
logger.info(f"Loaded samples from {config["long_read_sample_csv"]}")
long_read_samples = long_read_samples_df.set_index("sample").to_dict(orient="index")
logger.info(f"Analyzed {len(long_read_samples)} long-read samples from the CSV file.")