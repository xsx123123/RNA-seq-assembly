# --- 1. Load Libraries ---
library(clusterProfiler)
library(dplyr)
library(tidyr)      # Added for separate()
library(stringr)    # Added for str_replace()
library(log4r)
library(ggplot2)    # Added for ggplot()
library(forcats)    # Added for fct_reorder()

# --- 2. Logger Setup ---
# Create a logger object
logger <- create.logger()
# Set the log file path
logfile(logger) <- 'go_enrichment_script.log'
# Set the logging level (DEBUG will show all messages)
level(logger) <- 'DEBUG'

info(logger, "Script started. Logger initialized.")

# --- 3. Optimized Functions ---

#' Format GO Annotation Data
#' 
#' Efficiently splits the master GO data frame into TERM2GENE and TERM2NAME
#' lists for BP, MF, and CC ontologies.
#'
#' @param go_data A data frame with columns: GO_ID, GeneID, GO_Type, GO_Description
#' @return A list containing six data frames for clusterProfiler.
go_data_format <- function(go_data) {
  
  info(logger, "Starting GO data formatting.")
  
  # Optimization: Filter each ontology only ONCE
  
  # --- BP (Biological Process) ---
  info(logger, "Processing Biological Process (BP) data...")
  bp_data <- go_data |>
    dplyr::filter(GO_Type == "biological_process")
  
  bp_term2gene <- bp_data |>
    dplyr::select(GO_ID, GeneID) |>
    dplyr::distinct()
  bp_term2name <- bp_data |>
    dplyr::select(GO_ID, GO_Description) |>
    dplyr::distinct()
  debug(logger, paste("BP: Found", nrow(bp_term2gene), "TERM2GENE and", nrow(bp_term2name), "TERM2NAME mappings."))
  
  # --- MF (Molecular Function) ---
  info(logger, "Processing Molecular Function (MF) data...")
  mf_data <- go_data |>
    dplyr::filter(GO_Type == "molecular_function")
  
  mf_term2gene <- mf_data |>
    dplyr::select(GO_ID, GeneID) |>
    dplyr::distinct()
  mf_term2name <- mf_data |>
    dplyr::select(GO_ID, GO_Description) |>
    dplyr::distinct()
  debug(logger, paste("MF: Found", nrow(mf_term2gene), "TERM2GENE and", nrow(mf_term2name), "TERM2NAME mappings."))
  
  # --- CC (Cellular Component) ---
  info(logger, "Processing Cellular Component (CC) data...")
  cc_data <- go_data |>
    dplyr::filter(GO_Type == "cellular_component")
  
  cc_term2gene <- cc_data |>
    dplyr::select(GO_ID, GeneID) |>
    dplyr::distinct()
  cc_term2name <- cc_data |>
    dplyr::select(GO_ID, GO_Description) |>
    dplyr::distinct()
  debug(logger, paste("CC: Found", nrow(cc_term2gene), "TERM2GENE and", nrow(cc_term2name), "TERM2NAME mappings."))
  
  info(logger, "GO data formatting complete.")
  
  return(list(
    bp_term2gene = bp_term2gene,
    bp_term2name = bp_term2name,
    mf_term2gene = mf_term2gene,
    mf_term2name = mf_term2name,
    cc_term2gene = cc_term2gene,
    cc_term2name = cc_term2name
  ))
}


#' Format GO Descriptions for Plotting
#'
#' Vectorized function to wrap long GO term descriptions after the third word.
#'
#' @param description_vector A character vector of GO descriptions.
#' @return A character vector with formatted descriptions (containing '\n').
deal_description <- function(description_vector) {
  
  info(logger, "Formatting descriptions for plot labels (wrapping long text).")
  
  # Optimization: Replaced slow for-loop with a single, vectorized regex replacement
  new_descriptions <- stringr::str_replace(
    description_vector,
    # Regex: (3 words) + (space) + (rest of the string)
    pattern = "^(\\S+\\s+\\S+\\s+\\S+)\\s+(.*)$",
    # Replacement: (3 words) + (newline) + (rest of the string)
    replacement = "\\1\n\\2"
  )
  
  return(new_descriptions)
}


#' Run and Plot GO Enrichment
#'
#' Performs GO enrichment analysis for BP, MF, and CC, then generates
#' a dot plot of the top 10 results from each category (ranked by GeneRatio).
#'
#' @param go_dataset The formatted list object from go_data_format().
#' @param degs A character vector of differentially expressed gene IDs.
#' @return A list containing the ggplot object ($plot) and the three enrichResult objects.
GO_Enrichment <- function(go_dataset, degs) {
  
  info(logger, "Starting GO Enrichment process.")
  info(logger, paste("Input:", length(degs), "DEGs for enrichment."))
  
  # --- Helper Function (DRY Principle) ---
  # This function processes the enricher result, calculates numeric GeneRatio,
  # and filters for the Top 10 by GeneRatio.
  process_enrichment_result <- function(ego, ontology_name) {
    
    if (is.null(ego) || nrow(ego@result) == 0) {
      warn(logger, paste("No significant", ontology_name, "terms were found. Skipping."))
      return(NULL) # Return NULL if no results
    }
    
    debug(logger, paste("Processing", ontology_name, "results. Found", nrow(ego@result), "terms."))
    
    result_df <- ego |>
      as.data.frame() |>
      # BUG FIX: Convert to numeric *before* filtering/arranging
      tidyr::separate(col = GeneRatio, into = c("GR1", "GR2"), sep = "/", remove = FALSE) |>
      dplyr::mutate(GeneRatio_numeric = (as.numeric(GR1) / as.numeric(GR2))) |>
      # BUG FIX: Use slice_max to get top 10 by the numeric GeneRatio
      dplyr::slice_max(order_by = GeneRatio_numeric, n = 10) |>
      dplyr::mutate(ONTOLOGY = ontology_name)
    
    return(result_df)
  }
  
  # --- 1. Run Enrichment ---
  
  # BP
  info(logger, "Running enrichment for Biological Process (BP)...")
  bp_genes_in_go <- unique(go_dataset$bp_term2gene$GeneID)
  debug(logger, paste("BP Universe size:", length(bp_genes_in_go)))
  bp_ego <- enricher(
    gene = degs,
    TERM2GENE = go_dataset$bp_term2gene,
    TERM2NAME = go_dataset$bp_term2name,
    universe = bp_genes_in_go,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  
  # CC
  info(logger, "Running enrichment for Cellular Component (CC)...")
  cc_genes_in_go <- unique(go_dataset$cc_term2gene$GeneID)
  debug(logger, paste("CC Universe size:", length(cc_genes_in_go)))
  cc_ego <- enricher(
    gene = degs,
    TERM2GENE = go_dataset$cc_term2gene,
    TERM2NAME = go_dataset$cc_term2name,
    universe = cc_genes_in_go,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  
  # MF
  info(logger, "Running enrichment for Molecular Function (MF)...")
  mf_genes_in_go <- unique(go_dataset$mf_term2gene$GeneID)
  debug(logger, paste("MF Universe size:", length(mf_genes_in_go)))
  mf_ego <- enricher(
    gene = degs,
    TERM2GENE = go_dataset$mf_term2gene,
    TERM2NAME = go_dataset$mf_term2name,
    universe = mf_genes_in_go,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  
  info(logger, "Enrichment complete. Processing results for plotting.")
  
  # --- 2. Process and Bind Results ---
  
  bp_ego_result <- process_enrichment_result(bp_ego, "BP")
  mf_ego_result <- process_enrichment_result(mf_ego, "MF")
  cc_ego_result <- process_enrichment_result(cc_ego, "CC")
  
  # Bind all results together (NULLs will be automatically dropped)
  enrichment <- dplyr::bind_rows(bp_ego_result, mf_ego_result, cc_ego_result)
  
  # --- 3. Draw Plot ---
  
  if (nrow(enrichment) == 0) {
    error(logger, "No enrichment results found for any ontology. Cannot generate plot.")
    return(list(
      plot = NULL,
      bp_ego_res = bp_ego,
      cc_ego_res = cc_ego,
      mf_ego_res = mf_ego
    ))
  }
  
  info(logger, paste("Generating plot with", nrow(enrichment), "total terms."))
  
  # Apply description formatting
  enrichment$Description <- deal_description(enrichment$Description)
  
  # Reorder Description factor based on GeneRatio_numeric for plotting
  # This ensures the plot is sorted correctly
  p <- enrichment |>
    dplyr::mutate(Description = forcats::fct_reorder(Description, GeneRatio_numeric)) |>
    ggplot(aes(x = Description, y = GeneRatio_numeric)) + # Use the numeric value
    geom_point(aes(size = Count, color = qvalue)) +
    scale_color_gradient(low = "#f5af19", high = "#EF3B36") +
    coord_flip() +
    facet_grid(facets = "ONTOLOGY", scales = "free_y", space = "free_y") +
    labs(
      fill = "-log10 (Pvalue) ",
      x = "GO term",
      y = "GeneRatio", # Y-axis is now numeric
      title = "Enriched GO Terms"
    ) +
    theme_bw()
  
  info(logger, "GO Enrichment function finished.")
  
  return(list(
    plot = p,
    bp_ego_res = bp_ego,
    cc_ego_res = cc_ego,
    mf_ego_res = mf_ego
  ))
}

# --- 4. Example Usage (Mock Data) ---
# You would replace this with your actual data loading
#
# info(logger, "Creating mock data for demonstration.")
# 
# # Mock go_data
# go_data <- data.frame(
#   GO_ID = rep(c("GO:0006810", "GO:0007165", "GO:0005576", "GO:0003677"), each = 50),
#   GeneID = paste0("Gene", 1:200),
#   GO_Type = rep(c("biological_process", "biological_process", "cellular_component", "molecular_function"), each = 50),
#   GO_Description = rep(c("transporter activity very long description", 
#                          "signal transduction and more words", 
#                          "extracellular region or space", 
#                          "DNA binding and other things"), each = 50)
# )
# 
# # Mock DEGs
# degs <- paste0("Gene", 1:80)
# 
# info(logger, "Formatting mock GO data...")
# formatted_go <- go_data_format(go_data)
# 
# info(logger, "Running GO enrichment with mock data...")
# enrichment_results <- GO_Enrichment(go_dataset = formatted_go, degs = degs)
# 
# if (!is.null(enrichment_results$plot)) {
#   print(enrichment_results$plot)
#   info(logger, "Plot generated successfully.")
# } else {
#   warn(logger, "Plot could not be generated.")
# }
# 
# info(logger, "Script finished.")