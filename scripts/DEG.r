# --- 1. Load Libraries ---
# It's best practice to load all libraries at the top of your script.
library(log4r)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrepel)
library(tibble) # For rownames_to_column

# --- 2. Logger Setup ---
# Create a logger object
logger <- create.logger()
# Set the log file path
logfile(logger) <- 'deg_analysis.log'
# Set the logging level (DEBUG will show all messages)
level(logger) <- 'DEBUG'

info(logger, "--- Script started. Logger initialized. ---")

# --- 3. Optimized Functions ---

#' Check Metadata Structure
#'
#' Validates that the metadata file contains the required columns.
#'
#' @param metadata A data.frame for DESeq2 colData.
#' @param logger A log4r logger object.
Check_metadata <- function(metadata, logger) {
  
  info(logger, "Checking metadata integrity...")
  list <- colnames(metadata)
  
  # Optimization: Check for presence, not position (more robust)
  if (!'Sample' %in% list) {
    msg <- "Metadata check FAILED: 'Sample' column is missing. This column is required to match count data."
    error(logger, msg)
    stop(msg)
  }
  
  if (!'Group' %in% list) {
    msg <- "Metadata check FAILED: 'Group' column is missing. This column is required for the design formula."
    error(logger, msg)
    stop(msg)
  }
  
  # Optimization: Simplified factor check
  if (!is.factor(metadata$Group)) {
    warn(logger, "Metadata 'Group' column is NOT a factor. It will be converted to a factor, but setting levels manually in deg_by_deseq2 is recommended.")
  }
  
  info(logger, "Metadata check PASSED. 'Sample' and 'Group' columns are present.")
}


#' Run DEG Analysis using DESeq2
#'
#' Performs DESeq2 analysis, including data alignment and PCA plot generation.
#'
#' @param metadata The colData data.frame. Must contain 'Sample' and 'Group'.
#' @param data The raw count matrix. Column names must match 'Sample' IDs in metadata.
#' @param treat The name of the treatment group (numerator).
#' @param control The name of the control group (denominator).
#' @param pca_save_dir Directory to save PCA plots.
#' @param name Base name for output files.
#' @param logger A log4r logger object.
#' @return A list containing the DESeq results data.frame ($res) and the VST object ($vsd).
deg_by_deseq2 <- function(metadata,
                          data,
                          treat = 'CK_L_10',
                          control = 'CK_L_3',
                          pca_save_dir = "./",
                          name = 'CK_L',
                          logger) {
  
  info(logger, paste("Starting DESeq2 analysis for comparison:", name))
  debug(logger, paste("Treatment group:", treat, "| Control group:", control))
  
  # 1. Check metadata
  Check_metadata(metadata, logger)
  
  # 2. CRITICAL FIX: Align count data and metadata
  # This ensures colnames(data) and metadata$Sample are in the exact same order.
  info(logger, "Aligning count matrix columns to metadata 'Sample' order.")
  
  # Ensure 'Sample' is character for matching
  metadata$Sample <- as.character(metadata$Sample) 
  
  # Find common samples
  common_samples <- intersect(metadata$Sample, colnames(data))
  if (length(common_samples) == 0) {
    msg <- "No common samples found between metadata 'Sample' column and count data column names."
    error(logger, msg)
    stop(msg)
  }
  if (length(common_samples) < nrow(metadata)) {
    warn(logger, "Some samples in metadata were NOT found in count data.")
  }
  if (length(common_samples) < ncol(data)) {
    warn(logger, "Some samples in count data were NOT found in metadata.")
  }
  
  # Filter and reorder both
  metadata_aligned <- metadata |> 
    dplyr::filter(Sample %in% common_samples) |>
    dplyr::arrange(Sample)
  
  data_aligned <- data[, metadata_aligned$Sample]
  
  # Final check
  stopifnot(all(colnames(data_aligned) == metadata_aligned$Sample))
  debug(logger, paste("Alignment complete.", ncol(data_aligned), "samples will be used."))
  
  # 3. Create DESeqDataSet
  info(logger, "Creating DESeqDataSet object.")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(data_aligned), # Use aligned data
    colData = metadata_aligned,     # Use aligned metadata
    design = ~Group
  )
  
  # 4. Set factor levels for comparison
  info(logger, "Setting factor levels for DESeq comparison.")
  dds$Group <- factor(dds$Group, levels = c(control, treat))
  
  # 5. Pre-filtering
  dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]
  debug(logger, paste("Pre-filtering complete. Kept", nrow(dds), "genes."))
  
  # 6. Run DESeq
  info(logger, "Running DESeq(). This may take a moment...")
  dds <- DESeq2::DESeq(dds)
  
  # 7. Get results
  res <- DESeq2::results(dds, pAdjustMethod = "fdr", alpha = 0.05)
  res <- res[order(res$padj), ]
  info(logger, "DESeq analysis complete. Results extracted.")
  debug(logger, capture.output(summary(res)))
  
  # 8. PCA Plot
  info(logger, "Generating VST and PCA plot.")
  vsd <- DESeq2::vst(dds, blind = FALSE)
  pcaData <- DESeq2::plotPCA(vsd,
                             ntop = 500,
                             intgroup = c("Group"),
                             returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  p <- ggplot2::ggplot(pcaData, aes(PC1, PC2, color = Group)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_minimal() +
    guides(color = guide_legend(
      keywidth = 1, keyheight = 2, ncol = 1,
      override.aes = list(size = 3)
    ))
  
  # 9. Save PCA
  pca_png_path <- file.path(pca_save_dir, paste0(name, '-pca.png'))
  pca_pdf_path <- file.path(pca_save_dir, paste0(name, '-pca.pdf'))
  
  ggsave(pca_png_path, width = 5, height = 4, dpi = 300, plot = p)
  ggsave(pca_pdf_path, width = 5, height = 4, plot = p)
  info(logger, paste("PCA plots saved to:", pca_save_dir))
  
  # 10. Return results
  return(list(res = as.data.frame(res), vsd = vsd))
}


#' Draw a Volcano Plot
#'
#' Generates a publication-ready volcano plot from DESeq2 results.
#'
#' @param deg_result A data.frame from DESeq2 results. Must have gene IDs in rownames or a 'Symbol' column.
#' @param pvalCutoff Significance cutoff (e.g., 0.05).
#' @param LFCCutoff Log-fold-change cutoff (e.g., 1).
#' @param EXP_NAEE Experiment name string for the plot title.
#' @param y_aes Optional manual override for Y-axis limit.
#' @param p_value_type Column to use for p-value (default: "padj"). Can be "pvalue".
#' @param label_column Column to use for labels (default: "Symbol").
#' @param logger A log4r logger object.
#' @return A ggplot object.
DrawVolcano <- function(deg_result,
                        pvalCutoff,
                        LFCCutoff,
                        EXP_NAEE,
                        y_aes,
                        p_value_type = "padj",
                        label_column = "Symbol",
                        logger) {
  
  info(logger, paste("Drawing volcano plot for:", EXP_NAEE))
  debug(logger, paste("Settings: pvalCutoff =", pvalCutoff, "LFCCutoff =", LFCCutoff))
  debug(logger, paste("Using p-value column:", p_value_type, "| Label column:", label_column))
  
  # --- 1. Check for Label Column ---
  # If the specified label_column doesn't exist, use rownames.
  if (!label_column %in% colnames(deg_result)) {
    warn(logger, paste0("Label column '", label_column, "' not found. Using rownames for labels."))
    deg_result <- tibble::rownames_to_column(deg_result, var = "rowname_label")
    label_column <- "rowname_label"
  }
  
  # --- 2. Data Preparation ---
  
  # Check p_value_type exists
  if (!p_value_type %in% colnames(deg_result)) {
    msg <- paste("P-value column '", p_value_type, "' not found in deg_result. Available columns:", paste(colnames(deg_result), collapse=", "))
    error(logger, msg)
    stop(msg)
  }
  
  # Remove NA p-values
  deg_result <- deg_result |>
    dplyr::filter(!is.na(.data[[p_value_type]]))
  
  # Handle p-values of 0
  min_p_val <- .Machine$double.xmin
  if (min(deg_result[[p_value_type]], na.rm = TRUE) == 0) {
    warn(logger, "One or more p-values is 0. Converting to minimum possible value.")
    deg_result[[p_value_type]][which(deg_result[[p_value_type]] == 0)] <- min_p_val
  }
  
  # Create plotting columns
  deg_result <- deg_result |>
    dplyr::mutate(
      log10 = -log10(.data[[p_value_type]]),
      Group = "Non-significant"
    )
  
  # Add UP/DOWN/NO Symbol TAG
  deg_result$Group[which((deg_result[[p_value_type]] < pvalCutoff) & (deg_result$log2FoldChange > LFCCutoff))] <- "Up-regulated"
  deg_result$Group[which((deg_result[[p_value_type]] < pvalCutoff) & (deg_result$log2FoldChange < -LFCCutoff))] <- "Down-regulated"
  
  # --- 3. Prepare Labels (Top 15) ---
  
  # Sort by p-value to get top genes
  deg_result_sorted <- deg_result |> dplyr::arrange(.data[[p_value_type]])
  
  deg_result_up <- head(subset(deg_result_sorted, Group == "Up-regulated"), 15)
  deg_result_down <- head(subset(deg_result_sorted, Group == "Down-regulated"), 15)
  label_data <- dplyr::bind_rows(deg_result_up, deg_result_down)
  
  debug(logger, paste("Total significant: Up =", sum(deg_result$Group == "Up-regulated"), 
                      "| Down =", sum(deg_result$Group == "Down-regulated")))
  
  # --- 4. Dynamic Axis Limits (Your Logic) ---
  
  # Y-axis limit
  if (missing(y_aes)) {
    y_vals <- deg_result$log10[is.finite(deg_result$log10)]
    y_1 <- sort(y_vals, decreasing = TRUE)[1]
    y_2 <- sort(y_vals, decreasing = TRUE)[2]
    
    if (max(y_vals) > 300) {
      y_aes_value <- 250
    } else {
      if (!is.na(y_1) && !is.na(y_2) && y_1 / y_2 > 1.4) {
        y_aes_value <- (y_1 + y_2) / 2
      } else {
        y_aes_value <- max(y_vals) * 1.1
      }
    }
  } else {
    y_aes_value <- y_aes
    debug(logger, paste("Using manual y-axis limit:", y_aes_value))
  }
  
  # X-axis limit
  x_aes_vals <- deg_result$log2FoldChange
  x_aes_vals <- na.omit(x_aes_vals)
  if (max(x_aes_vals) > 7.5) {
    x_aes_limit <- 7.5
  } else {
    x_aes_limit <- max(x_aes_vals)
  }
  x_aes_limit <- x_aes_limit * 1.2
  
  # --- 5. Draw Volcano Plot ---
  
  # OPTIMIZATION: Simplified plot using aes(color=Group) instead of 5 geoms
  p <- ggplot2::ggplot(deg_result, aes(x = log2FoldChange, y = log10)) +
    
    # Single geom_point layer for all points
    geom_point(aes(color = Group), size = 0.5, alpha = 0.6) +
    
    # Manual color scale
    scale_color_manual(values = c(
      "Up-regulated" = "#e41749",
      "Down-regulated" = "#41b6e6",
      "Non-significant" = "#C7C7C7"
    )) +
    
    # Threshold lines
    # OPTIMIZATION: Dynamic y-intercept based on pvalCutoff
    geom_vline(xintercept = LFCCutoff, lty = 2, col = "black", lwd = 0.1) +
    geom_vline(xintercept = -LFCCutoff, lty = 2, col = "black", lwd = 0.1) +
    geom_hline(yintercept = -log10(pvalCutoff), lty = 2, col = "black", lwd = 0.1) +
    
    # Labels
    labs(
      x = bquote("RNA-seq " * log[2] * " fold change " * .(EXP_NAEE) * ""),
      y = bquote(-log[10] ~ .(p_value_type)), # Dynamic y-axis label
      title = paste0(EXP_NAEE, " Volcano Plot")
    ) +
    
    # Text labels for top genes
    ggrepel::geom_text_repel(
      data = label_data,
      aes(label = .data[[label_column]]), # Use .data[[]]
      size = 0.7, colour = "black", fontface = "bold.italic",
      segment.alpha = 0.5, segment.size = 0.15, segment.color = "black",
      min.segment.length = 0, box.padding = unit(0.2, "lines"),
      point.padding = unit(0, "lines"), force = 20, max.iter = 3e3,
      max.overlaps = 25, arrow = arrow(length = unit(0.02, "inches"))
    ) +
    
    # Axis scales
    scale_x_continuous(limits = c(-x_aes_limit, x_aes_limit), n.breaks = 8) +
    scale_y_continuous(limits = c(0, y_aes_value), n.breaks = 10) +
    
    # Theme (Kept your custom theme)
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 5, family = "sans", face = "bold"),
      legend.position = "none",
      text = element_text(size = 4, family = "sans"),
      title = element_text(size = 4),
      axis.title.x = element_text(size = 3, family = "sans", face = "bold"),
      axis.title.y = element_text(size = 3, family = "sans", face = "bold"),
      axis.text.x = element_text(size = 2.5, family = "sans"),
      axis.text.y = element_text(size = 2.5, family = "sans"),
      axis.line.x = element_line(linetype = 1, color = "#606c70", size = 0.2),
      axis.line.y = element_line(linetype = 1, color = "#606c70", size = 0.2),
      axis.ticks.x = element_line(color = "#606c70", size = 0.15, lineend = 0.05),
      axis.ticks.length = unit(.08, "lines"),
      axis.ticks.y = element_line(color = "#606c70", size = 0.15, lineend = 0.05)
    )
  
  info(logger, "Volcano plot generated successfully.")
  return(p)
}

# --- 5. Example Usage (Mock Data) ---
#
# info(logger, "Creating mock data for demonstration.")
# 
# # Mock count data
# set.seed(123)
# counts <- matrix(rpois(1000 * 6, lambda = 10), ncol = 6)
# colnames(counts) <- c("ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3")
# rownames(counts) <- paste0("Gene", 1:1000)
# 
# # Add some differential expression
# counts[1:50, 4:6] <- counts[1:50, 4:6] * 10 # Up-regulated
# counts[51:100, 4:6] <- counts[51:100, 4:6] / 10 # Down-regulated
# 
# # Mock metadata
# metadata <- data.frame(
#   Sample = c("ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"),
#   Group = factor(c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"))
# )
# 
# info(logger, "Running optimized deg_by_deseq2 function...")
# 
# # Run DESeq2
# deg_results_list <- deg_by_deseq2(
#   metadata = metadata,
#   data = counts,
#   treat = "Treatment",
#   control = "Control",
#   name = "Test_Comparison",
#   logger = logger
# )
# 
# # Add a mock 'Symbol' column for plotting
# deg_results_with_symbols <- deg_results_list$res
# deg_results_with_symbols$Symbol <- rownames(deg_results_with_symbols)
# 
# info(logger, "Running optimized DrawVolcano function...")
# 
# # Run Volcano Plot
# volcano_plot <- DrawVolcano(
#   deg_result = deg_results_with_symbols,
#   pvalCutoff = 0.05,
#   LFCCutoff = 1,
#   EXP_NAEE = "Treatment vs Control",
#   p_value_type = "padj",      # Use adjusted p-value
#   label_column = "Symbol",    # Use the column we just added
#   logger = logger
# )
# 
# # print(volcano_plot)
# 
# # ggsave("test_volcano.png", plot = volcano_plot, width = 6, height = 5, dpi = 300)
# info(logger, "--- Script finished. ---")