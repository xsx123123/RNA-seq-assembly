import pandas as pd
import numpy as np
import argparse
import sys
from loguru import logger
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Import plotnine library
from plotnine import (
    ggplot,
    aes,
    geom_line,
    geom_point,
    geom_vline,
    geom_area,
    labs,
    theme_bw,
    theme,
    element_text,
    element_blank,
    scale_x_continuous,
    scale_y_continuous,
    annotate,
    save_as_pdf_pages,
)

# Configure loguru logger
logger.remove()
logger.add(sys.stderr, format="<green>{time:HH:mm:ss}</green> | <level>{level}</level> | {message}", level="INFO")

def create_exn50_plot(input_file: str, output_prefix: str, formats: list, width: float, height: float):
    """
    Generates ExN50 plots and saves them separately or combined.

    Args:
        input_file: Path to the ExN50 statistics file.
        output_prefix: Prefix name for the output files.
        formats: List of output file formats (e.g., ['png', 'pdf']).
        width: Plot width (inches).
        height: Total plot height (inches).
    """
    logger.info(f"Starting to process file: {input_file}")

    # --- 1. Data Reading and Preparation ---
    try:
        exn50_data = pd.read_csv(input_file, sep='\t')
        
        required_cols = ['Ex', 'ExN50', 'num_transcripts']
        if not all(col in exn50_data.columns for col in required_cols):
            logger.error(f"Input file is missing required columns. File must contain: {required_cols}")
            sys.exit(1)

    except FileNotFoundError:
        logger.error(f"Error: Input file not found at {input_file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error occurred while reading  {e}")
        sys.exit(1)

    logger.info(f"Data successfully loaded. First 5 rows:\n{exn50_data.head().to_string()}")
    
    # Extract key values: E90N50 and E90 transcript count
    e90_row = exn50_data[exn50_data['Ex'] == 90]
    
    if e90_row.empty:
        logger.warning("No data row found where Ex=90. Skipping E90N50 annotation.")
        e90n50_value = None
        e90_transcripts = None
    else:
        e90n50_value = e90_row['ExN50'].iloc[0]
        e90_transcripts = e90_row['num_transcripts'].iloc[0]
        logger.info(f"Extracted key metrics: E90N50 = {e90n50_value} bp, E90 transcripts = {e90_transcripts:,}")

    # --- 2. Plotting ExN50 Curve (p1) ---
    max_exn50 = exn50_data['ExN50'].max()
    p1 = (
        ggplot(exn50_data, aes(x='Ex', y='ExN50'))
        + geom_line(color="#2E86AB", size=1.5)
        + geom_point(color="#2E86AB", size=2, alpha=0.6)
        
        + geom_vline(xintercept=90, linetype="dashed", color="#A23B72", size=1)
        
        # Add annotation for E90N50 value (if present)
        + (
            annotate("text",
                     x=90,
                     y=max_exn50 * 1.05,
                     label=f"E90N50 = {e90n50_value} bp",
                     color="#A23B72",
                     size=11, 
                     ha='center',
                     fontweight='bold')
            if e90n50_value is not None else []
        )
        
        + scale_x_continuous(breaks=np.arange(0, 101, 10))
        + labs(
            title="ExN50 Statistics for Transcriptome Assembly",
            x="Top expressed transcripts (%)",
            y="ExN50 contig length (bp)"
        )
        + theme_bw()
        + theme(
            plot_title=element_text(hjust=0.5, face="bold", size=14),
            axis_title=element_text(size=12, face="bold"),
            axis_text=element_text(size=11),
            panel_grid_minor=element_blank(),
            figure_size=(width, height / 2.2)
        )
    )

    # --- 3. Plotting Transcript Count Curve (p2) ---
    max_num_transcripts = exn50_data['num_transcripts'].max()
    p2 = (
        ggplot(exn50_data, aes(x='Ex', y='num_transcripts'))
        + geom_area(fill="#F18F01", alpha=0.4)
        + geom_line(color="#C73E1D", size=1.5)
        + geom_point(color="#C73E1D", size=2, alpha=0.6)
        
        + geom_vline(xintercept=90, linetype="dashed", color="#A23B72", size=1)

        # Add annotation for E90 transcript count (if present)
        + (
            annotate("text",
                     x=90,
                     y=max_num_transcripts * 1.05,
                     label=f"E90 transcripts = {e90_transcripts:,}",
                     color="#A23B72",
                     size=11,
                     ha='center',
                     fontweight='bold')
            if e90_transcripts is not None else []
        )
        
        + scale_x_continuous(breaks=np.arange(0, 101, 10))
        # Use lambda function with f-string formatting for comma separation
        + scale_y_continuous(labels=lambda l: [f"{int(x):,}" for x in l])
        + labs(
            title="Number of Transcripts by Expression Percentile",
            x="Top expressed transcripts (%)",
            y="Number of transcripts"
        )
        + theme_bw()
        + theme(
            plot_title=element_text(hjust=0.5, face="bold", size=14),
            axis_title=element_text(size=12, face="bold"),
            axis_text=element_text(size=11),
            panel_grid_minor=element_blank(),
            figure_size=(width, height / 2.2)
        )
    )

    # --- 4. 保存单独的图表，然后用 matplotlib 组合 ---
    import tempfile
    import os
    from PIL import Image
    
    # 创建临时文件保存单独的图
    with tempfile.TemporaryDirectory() as tmpdir:
        temp_p1 = os.path.join(tmpdir, "p1_temp.png")
        temp_p2 = os.path.join(tmpdir, "p2_temp.png")
        
        # 保存两个图到临时文件
        p1.save(temp_p1, dpi=300, verbose=False)
        p2.save(temp_p2, dpi=300, verbose=False)
        
        # 使用 matplotlib 组合图片
        fig, axes = plt.subplots(2, 1, figsize=(width, height))
        
        # 读取并显示图片
        img1 = Image.open(temp_p1)
        img2 = Image.open(temp_p2)
        
        axes[0].imshow(img1)
        axes[0].axis('off')
        
        axes[1].imshow(img2)
        axes[1].axis('off')
        
        # 添加总标题
        fig.suptitle(
            'Transcriptome Assembly Quality Assessment',
            fontsize=16,
            fontweight='bold',
            y=0.995
        )
        
        # 调整布局
        plt.tight_layout(rect=[0, 0, 1, 0.99])
        
        # --- 5. 保存组合图 ---
        for fmt in formats:
            filename = f"{output_prefix}.{fmt.lower()}"
            try:
                fig.savefig(filename, dpi=300, bbox_inches='tight')
                logger.info(f"Combined plot successfully saved to: {filename}")
            except Exception as e:
                logger.error(f"Error occurred while saving file {filename}: {e}")
        
        plt.close(fig)
    
    # 同时保存单独的图
    p1.save(f"{output_prefix}_ExN50.png", dpi=300, verbose=False)
    p2.save(f"{output_prefix}_transcripts.png", dpi=300, verbose=False)
    logger.info(f"Individual plots saved: {output_prefix}_ExN50.png, {output_prefix}_transcripts.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plots ExN50 statistics curves for transcriptome assembly.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        '-i', '--input', 
        required=True, 
        help="Path to the ExN50 statistics file (e.g., ExN50.transcript.stats)."
    )
    parser.add_argument(
        '-o', '--output_prefix', 
        default='ExN50_combined',
        help="Prefix name for the output files (default: ExN50_combined)."
    )
    parser.add_argument(
        '-F', '--formats', 
        nargs='+', 
        default=['pdf', 'png'],
        choices=['pdf', 'png', 'svg', 'jpg', 'jpeg'],
        help="List of output file formats (space-separated, e.g., -F png pdf)."
    )
    parser.add_argument(
        '--width', 
        type=float, 
        default=10, 
        help="Plot width (inches) (default: 10)."
    )
    parser.add_argument(
        '--height', 
        type=float, 
        default=12, 
        help="Plot height (inches) (default: 12)."
    )

    parser.add_argument(
        '--simulate',
        action='store_true',
        help=argparse.SUPPRESS 
    )

    args = parser.parse_args()

    if args.simulate:
        logger.warning("Using simulation mode to create test file 'dummy_exn50.stats'.")
        input_file = "dummy_exn50.stats"
        dummy_data = {
            'Ex': np.arange(1, 101, 1),
            'ExN50': np.random.randint(1500, 2500, 100),
            'num_transcripts': np.random.randint(100, 200000, 100)
        }
        pd.DataFrame(dummy_data).to_csv(input_file, sep='\t', index=False)
        logger.info(f"Test file '{input_file}' created.")
        args.input = input_file
    
    create_exn50_plot(
        args.input, 
        args.output_prefix, 
        args.formats, 
        args.width, 
        args.height
    )
