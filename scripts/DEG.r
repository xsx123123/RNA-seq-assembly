#!/usr/bin/env Rscript

# ================================================================= #
# Script: run_deg_analysis.R
# Author: zhang jian (refactored by Hakimi)
# Date: 2025-10-02
# Version: 2.0v
# Description: A general-purpose script for DEG analysis from 
#              Salmon quantification results using DESeq2. 
#              Accepts command-line arguments.
# ================================================================= #

# --- 1. Load Packages --- #
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrepel"))

# ================================================================= #
# --- 2. Define Functions (Your original functions) --- #
# ================================================================= #

Check_metadata <- function(metadata){
  list <- colnames(metadata)
  if ('Sample' %in% list){
    cat('  - Metadata check: `Sample` column found. OK.\n')
  } else {
    stop('Metadata check failed: `Sample` column is missing.')
  }
  if ('Group' %in% list){
    cat('  - Metadata check: `Group` column found. OK.\n')
  } else {
    stop('Metadata check failed: `Group` column is missing.')
  }
  if (!is.factor(metadata$Group)){
      cat('  - Metadata check: `Group` column is not a factor. Converting...\n')
      metadata$Group <- as.factor(metadata$Group)
  } else {
      cat('  - Metadata check: `Group` column is already a factor. OK.\n')
  }
  return(metadata)
}

deg_by_deseq2 <- function(metadata,
                          data,
                          treat,
                          control,
                          pca_save_dir = "./",
                          name = 'comparison'
){
  # check metadata
  metadata <- Check_metadata(metadata)
  
  # create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = round(data), 
                                colData = metadata, 
                                design = ~Group)
  
  # set factor levels for comparison
  dds$Group <- factor(dds$Group, levels = c(control, treat))
  
  # pre-filtering
  dds <- dds[rowSums(counts(dds)) > 1, ]
  
  # run DESeq2
  cat(paste0("  - Running DESeq2 for: ", name, " (", treat, " vs ", control, ")\n"))
  dds <- DESeq(dds)
  res <- results(dds, pAdjustMethod = "fdr", alpha = 0.05)
  res <- res[order(res$padj),]
  summary(res)
  
  # --- PCA Plot --- #
  vsd <- vst(dds, blind = FALSE)
  pcaData <- plotPCA(vsd, 
                     ntop = 500,
                     intgroup = c("Group"),
                     returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  p <- ggplot(pcaData, aes(PC1, PC2, color = Group)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title = paste("PCA Plot for", name)) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5))

  # Save PCA plot
  ggsave(file.path(pca_save_dir, paste0(name, '-pca.png')),
         width = 5, height = 4, dpi = 300, plot = p)
  ggsave(file.path(pca_save_dir, paste0(name, '-pca.pdf')),
         width = 5, height = 4, plot = p)
  
  cat(paste0("  - PCA plot saved to: ", pca_save_dir, "\n"))
  
  return(res)
}

DrawVolcano <- function(deg_result, pvalCutoff, LFCCutoff, EXP_NAME, y_aes_max){
  # Rename pvalue to avoid conflict with function name
  deg_result <- deg_result %>% dplyr::rename(p_val_adj = padj)
  
  # Handle cases where p-values might be NA
  deg_result <- deg_result[!is.na(deg_result$p_val_adj), ]
  
  # Create a column for -log10(p-adj)
  deg_result <- deg_result %>% dplyr::mutate(log10_padj = -log10(p_val_adj))

  # Add UP/DOWN/NO Symbol TAG
  deg_result$Group <- "Non-significant"
  deg_result$Group[which((deg_result$p_val_adj < pvalCutoff) & (deg_result$log2FoldChange > LFCCutoff))] <- "Up-regulated"
  deg_result$Group[which((deg_result$p_val_adj < pvalCutoff) & (deg_result$log2FoldChange < -LFCCutoff))] <- "Down-regulated"
  
  # Sort by adjusted p-value to label top genes
  deg_result <- deg_result %>% dplyr::arrange(p_val_adj)
  
  # Get top 15 up and down genes for labeling
  deg_result_up <- head(subset(deg_result, Group == "Up-regulated"), 15)
  deg_result_down <- head(subset(deg_result, Group == "Down-regulated"), 15)
  
  # --- Axis limits --- #
  y_limit <- if (!missing(y_aes_max)) {
      y_aes_max
  } else {
      max_y <- max(deg_result$log10_padj[is.finite(deg_result$log10_padj)], na.rm = TRUE)
      if (max_y > 200) 200 else ceiling(max_y * 1.1)
  }
  
  x_limit_val <- max(abs(deg_result$log2FoldChange), na.rm = TRUE)
  x_limit <- if (x_limit_val > 10) 10 else ceiling(x_limit_val * 1.1)

  # --- Draw Volcano Plot --- #
  p <- ggplot(deg_result, aes(x = log2FoldChange, y = log10_padj)) +
    geom_point(aes(color = Group), size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("Down-regulated" = "#41b6e6", "Non-significant" = "#C7C7C7", "Up-regulated" = "#e41749")) +
    geom_vline(xintercept = c(-LFCCutoff, LFCCutoff), lty = 2, col = "black", lwd = 0.5) +
    geom_hline(yintercept = -log10(pvalCutoff), lty = 2, col = "black", lwd = 0.5) +
    labs(
        x = bquote("log"[2]*" Fold Change"),
        y = bquote("-log"[10]*" (Adjusted P-value)"),
        title = paste(EXP_NAME, "Volcano Plot")
    ) +
    geom_text_repel(
        data = rbind(deg_result_up, deg_result_down),
        aes(label = Symbol),
        size = 2.5,
        box.padding = 0.3,
        point.padding = 0.3,
        segment.alpha = 0.6,
        max.overlaps = 20
    ) +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    scale_y_continuous(limits = c(0, y_limit)) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank()
    )
    
  return(p)
}


# ================================================================= #
# --- 3. Command-Line Argument Parsing --- #
# ================================================================= #

parser <- ArgumentParser(description = 'Run DEG analysis using DESeq2.')
parser$add_argument('--salmon_dir', type = 'character', required = TRUE, help = 'Path to the directory containing Salmon quantification folders.')
parser$add_argument('--sample_info', type = 'character', required = TRUE, help = 'Path to the sample metadata CSV file. Must contain "Sample", "Group", and "ComparisonPair" columns.')
parser$add_argument('--transcript_ids', type = 'character', required = TRUE, help = 'Path to the text file with transcript IDs to keep.')
parser$add_argument('--output_dir', type = 'character', required = TRUE, help = 'Path to the directory where results will be saved.')
parser$add_argument('--pval', type = 'double', default = 0.05, help = 'Adjusted p-value cutoff for volcano plot.')
parser$add_argument('--lfc', type = 'double', default = 1.0, help = 'Log2 fold change cutoff for volcano plot. Default is 1 (meaning a 2-fold change).')

args <- parser$parse_args()

# Create output directory if it doesn't exist
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
}

# ================================================================= #
# --- 4. Main Workflow --- #
# ================================================================= #

cat("Starting DEG Analysis Workflow...\n")
cat("---------------------------------\n")

# --- 4.1. Load and Prepare Data --- #
cat("1. Loading data with tximport...\n")
sample_info <- fread(args$sample_info)
quant_files <- file.path(args$salmon_dir, sample_info$Sample, "quant.sf")
names(quant_files) <- sample_info$Sample

# Check if all files exist
if (!all(file.exists(quant_files))) {
  stop("Some quant.sf files are missing. Please check your --salmon_dir and --sample_info file.")
}

txi <- tximport(quant_files, type = "salmon", txOut = TRUE)
counts_matrix <- txi$counts

cat("2. Filtering for E90 transcriptome IDs...\n")
transcript_ids_to_keep <- fread(args$transcript_ids, header = FALSE)$V1
filtered_matrix <- as.data.frame(counts_matrix) %>%
  rownames_to_column(var = 'trans_name') %>%
  filter(trans_name %in% transcript_ids_to_keep) %>%
  column_to_rownames(var = 'trans_name')

cat(paste0("  - Original transcripts: ", nrow(counts_matrix), "\n"))
cat(paste0("  - Filtered transcripts: ", nrow(filtered_matrix), "\n"))

# --- 4.2. Loop Through Comparisons --- #
cat("3. Starting DEG analysis for each comparison pair...\n")

# Get unique comparison pairs from the sample info file
comparison_pairs <- unique(sample_info$ComparisonPair)

for (pair in comparison_pairs) {
  
  cat(paste0("\nProcessing comparison: ", pair, "\n"))
  
  # Subset metadata for the current comparison
  current_metadata <- sample_info %>% filter(ComparisonPair == pair) %>% as.data.frame()
  rownames(current_metadata) <- current_metadata$Sample
  
  # Subset the counts matrix
  current_counts <- filtered_matrix[, current_metadata$Sample]
  
  # Determine treat and control groups. We assume the group with a higher number is the treatment.
  # This logic can be adjusted if needed. For "CK_L_10" vs "CK_L_3", "CK_L_10" is treat.
  groups <- unique(current_metadata$Group)
  group_numbers <- as.numeric(gsub(".*_(\\d+)$", "\\1", groups))
  
  if(length(groups) != 2 || any(is.na(group_numbers))) {
      warning(paste("Could not automatically determine treat/control for", pair, ". Skipping."))
      next
  }

  treat_group <- groups[which.max(group_numbers)]
  control_group <- groups[which.min(group_numbers)]
  
  # Run DESeq2 analysis using your function
  deg_result <- deg_by_deseq2(
    metadata = current_metadata,
    data = current_counts,
    treat = treat_group,
    control = control_group,
    pca_save_dir = args$output_dir,
    name = pair
  )
  
  # Save DEG results to CSV
  deg_df <- as.data.frame(deg_result) %>% rownames_to_column(var = 'Symbol')
  output_csv_path <- file.path(args$output_dir, paste0(pair, "_deg_results.csv"))
  write.csv(deg_df, output_csv_path, row.names = FALSE)
  cat(paste0("  - DEG results saved to: ", output_csv_path, "\n"))
  
  # Draw and save Volcano plot
  volcano_plot <- DrawVolcano(
    deg_result = deg_df,
    pvalCutoff = args$pval,
    LFCCutoff = args$lfc,
    EXP_NAME = pair
  )
  
  output_volcano_png <- file.path(args$output_dir, paste0(pair, "_volcano.png"))
  output_volcano_pdf <- file.path(args$output_dir, paste0(pair, "_volcano.pdf"))
  
  ggsave(output_volcano_png, width = 7, height = 6, dpi = 300, plot = volcano_plot)
  ggsave(output_volcano_pdf, width = 7, height = 6, plot = volcano_plot)
  cat(paste0("  - Volcano plot saved to: ", args$output_dir, "\n"))
}

cat("\n---------------------------------\n")
cat("Workflow finished successfully! 🎉\n")