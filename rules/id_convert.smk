#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
from rich import inspect,print
# Load the samples CSV file into a DataFrame
samples_df = pd.read_csv(config["sample_csv"])
logger.info(f"Loaded samples from {config["sample_csv"]}")
# Convert the DataFrame to a dictionary for easy access
load_samples = samples_df.set_index("sample").to_dict(orient="index")
logger.info(f"Analyzed {len(load_samples)} samples from the CSV file.")