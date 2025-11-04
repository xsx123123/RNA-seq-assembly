#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import argparse
import hashlib
import logging
from typing import List
import concurrent.futures
from datetime import datetime

# --- Setup Logging ---
def setup_logger() -> logging.Logger:
    """Configures and returns a logger instance that writes to both console and a file."""
    logger = logging.getLogger("MD5_Checker")
    logger.setLevel(logging.INFO)  # Set the minimum logging level to INFO

    # Avoid adding handlers multiple times if the logger is already configured
    if not logger.handlers:
        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        
        # 1. Console Handler
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

        # 2. File Handler (writes logs to md5_checker.log)
        file_handler = logging.FileHandler("md5_checker.log", mode='w', encoding='utf-8')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        
    return logger

# Get the configured logger instance
logger = setup_logger()


# --- Core Functions ---
def calculate_md5(filepath: str, chunk_size: int = 8192) -> str:
    """
    Calculates the MD5 hash of a given file.
    The file is read in chunks to handle large files without consuming too much memory.

    Args:
        filepath (str): The path to the file.
        chunk_size (int): The size of each chunk to read in bytes.

    Returns:
        str: The lowercase hexadecimal MD5 hash of the file.
    """
    md5_hash = hashlib.md5()
    try:
        with open(filepath, "rb") as f:
            while chunk := f.read(chunk_size):
                md5_hash.update(chunk)
    except IOError as e:
        logger.error(f"Could not read file {filepath}: {e}")
        return ""
    return md5_hash.hexdigest()


def verify_file_task(task_info: tuple) -> dict:
    """
    Worker thread function to verify a single file's MD5 and return a result dictionary.
    
    Args:
        task_info (tuple): A tuple containing (expected_md5, file_to_check, relative_filename).
    
    Returns:
        dict: A dictionary containing detailed verification results.
    """
    expected_md5, file_to_check, _ = task_info
    result = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "full_path": file_to_check,
        "expected_md5": expected_md5,
        "actual_md5": "N/A",
        "status": "FAIL",
        "message": ""
    }

    if not os.path.exists(file_to_check):
        result["message"] = "File not found"
        return result

    actual_md5 = calculate_md5(file_to_check)
    if not actual_md5:
        result["message"] = "Could not calculate MD5 (read error)"
        return result
    
    result["actual_md5"] = actual_md5

    if actual_md5 == expected_md5:
        result["status"] = "PASS"
        result["message"] = "MD5 match"
    else:
        result["message"] = "MD5 mismatch"

    return result

def generate_report(results: List[dict], output_file: str):
    """Generates a TSV report from the verification results."""
    logger.info(f"--- Generating verification report to: {output_file} ---")
    try:
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            # Write header
            header = "CheckTime\tFilePath\tExpectedMD5\tActualMD5\tStatus\tMessage\n"
            f.write(header)
            # Write data rows
            for res in results:
                line = (
                    f"{res['timestamp']}\t"
                    f"{res['full_path']}\t"
                    f"{res['expected_md5']}\t"
                    f"{res['actual_md5']}\t"
                    f"{res['status']}\t"
                    f"{res['message']}\n"
                )
                f.write(line)
        logger.info("Report generation complete.")
    except IOError as e:
        logger.error(f"Could not write to report file {output_file}: {e}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Automatically finds MD5 checksum files in specified directories and verifies file integrity using multiple threads.",
        epilog="Example: python md5_checker.py /path/to/data1 /path/to/data2 -f checksums.txt -t 16 -o report.tsv"
    )
    parser.add_argument(
        "directories",
        nargs='+',
        help="One or more directory paths containing data files and an MD5 checksum file."
    )
    parser.add_argument(
        "-f", "--filename",
        type=str,
        default="MD5.txt",
        help="The name of the MD5 checksum file (default: MD5.txt)."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of threads to use for concurrent verification (default: 8)."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        help="Path to generate a TSV-formatted report file."
    )
    args = parser.parse_args()

    # --- Step 1: Discover and parse all MD5 files to create a task list ---
    tasks = []
    logger.info(f"--- Start searching for {args.filename} files and building task list ---")
    for directory in args.directories:
        if not os.path.isdir(directory):
            logger.warning(f"Provided path is not a directory, skipping: {directory}")
            continue
        
        md5_file_path = os.path.join(directory, args.filename)
        if not os.path.exists(md5_file_path):
            logger.warning(f"Could not find {args.filename} in directory {directory}, skipping.")
            continue
        
        logger.info(f"Found and parsing: {md5_file_path}")
        with open(md5_file_path, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f, 1):
                line = line.strip()
                if not line: continue

                parts = line.split(maxsplit=1)
                if len(parts) != 2:
                    logger.warning(f"Warning: Malformed line #{i} in {md5_file_path}, skipping: '{line}'")
                    continue
                
                expected_md5, relative_filename = parts
                file_to_check = os.path.join(directory, relative_filename)
                tasks.append((expected_md5, file_to_check, relative_filename))

    # --- Step 2: Execute verification tasks using a thread pool ---
    num_tasks = len(tasks)
    if num_tasks == 0:
        logger.warning("No files to verify were found. Script finished.")
        sys.exit(0)

    logger.info(f"--- Found {num_tasks} files to verify, starting {args.threads} threads for processing ---")
    
    results = []
    has_failures = False
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_task = {executor.submit(verify_file_task, task): task for task in tasks}

        for future in concurrent.futures.as_completed(future_to_task):
            relative_filename = future_to_task[future][2]
            try:
                result = future.result()
                results.append(result)
                
                if result['status'] == 'PASS':
                    logger.info(f"  [SUCCESS] MD5 match: {result['expected_md5']}  {relative_filename}")
                else:
                    has_failures = True
                    logger.error(f"  [FAIL] {result['message']}: {relative_filename}")
                    if result['message'] == 'MD5 mismatch':
                        logger.error(f"      Expected: {result['expected_md5']}")
                        logger.error(f"      Actual:   {result['actual_md5']}")
            except Exception as exc:
                has_failures = True
                logger.critical(f"An unexpected error occurred while processing {relative_filename}: {exc}")

    # --- Step 3: Generate report and output final results ---
    if args.output:
        generate_report(results, args.output)

    logger.info("======================================================")
    if has_failures:
        logger.critical("Errors were found during verification. Please check logs and the report file for details.")
        logger.info("======================================================")
        sys.exit(1)
    else:
        successful_checks = len(results)
        logger.info(f"All {successful_checks} files passed verification successfully!")
        logger.info("======================================================")

# --- Main script entry point ---
if __name__ == "__main__":
    main()