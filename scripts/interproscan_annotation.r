# --- 1. Load required R packages ---
# On first run, you may need to uncomment the following lines to install packages
# install.packages(c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "forcats", "UpSetR", "patchwork"))
# 
# # GO.db is a Bioconductor package, installed differently:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(c("GO.db", "AnnotationDbi"))

# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)     # For text processing
library(forcats)     # For reordering in ggplot2
library(UpSetR)      # For Upset plot
library(GO.db)       # For GO ontology
library(AnnotationDbi) # For GO ontology
library(patchwork)   # --- NEW: For combining ggplot plots ---

# Define input file
input_file <- "/data/jzhang/project/RNA-seq-assembly_2025.8.25/05.transcript_annotation/TD2_pep_interproscan_annotation/rnabloom_transcript_LongOrfs_ann.summary"
# Define output directory
output_dir <- "/data/jzhang/project/RNA-seq-assembly_2025.8.25/05.transcript_annotation/TD2_pep_interproscan_annotation"

# --- 2. Create output directory ---
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
cat(paste("Visualization results will be saved to:", output_dir, "\n"))

# --- 3. Load and preprocess data ---
cat("Loading data...\n")
# Read all columns as character to prevent issues with '-'
data <- read_tsv(input_file, col_types = cols(.default = "c"))

# Create a statistics data frame for Upset plot and pie chart
data_stats <- data %>%
  mutate(
    has_Pfam = ifelse(Pfam_Annotation != "-", 1, 0),
    has_PANTHER = ifelse(PANTHER_Annotation != "-", 1, 0),
    has_InterPro = ifelse(InterPro_Entries != "-", 1, 0),
    has_GO = ifelse(GO_Terms != "-", 1, 0)
  )

# --- 4. Generate plots ---

# === Plot 1: Overall Annotation Rate Pie Chart ===
cat("Generating... 1. Overall Annotation Rate Pie Chart\n")
total_proteins <- nrow(data_stats)
annotated_count <- data_stats %>%
  filter(has_Pfam == 1 | has_PANTHER == 1 | has_InterPro == 1 | has_GO == 1) %>%
  nrow()

stats_df <- data.frame(
  Category = c("Annotated", "Unannotated"),
  Count = c(annotated_count, total_proteins - annotated_count)
) %>%
  mutate(Percentage = Count / sum(Count),
         Label = paste0(Category, "\n(", scales::percent(Percentage, accuracy = 0.1), ")"))

p1 <- ggplot(stats_df, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), color = "white", fontface = "bold") +
  labs(title = "Overall Annotation Status",
       subtitle = paste("Total Proteins:", total_proteins)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12))

# --- Modify: Remove ggsave, will combine later ---
# ggsave(file.path(output_dir, "1_annotation_pie_chart.pdf"), p1, width = 8, height = 6)


# === Plot 2: Annotation Source Upset Plot ===
cat("Generating... 2. Annotation Source Upset Plot (saving separately)\n")
# UpsetR requires a data.frame
upset_df <- as.data.frame(data_stats[, c("has_Pfam", "has_PANTHER", "has_InterPro", "has_GO")])
# Rename for plotting
colnames(upset_df) <- c("Pfam", "PANTHER", "InterPro", "GO")

# Save as PDF (UpsetR cannot be combined with patchwork, so save separately)
pdf(file.path(output_dir, "2_annotation_upset_plot.pdf"), width = 10, height = 7)
upset(upset_df,
      sets = c("GO", "InterPro", "Pfam", "PANTHER"), # Reorder sets
      nintersects = 20, # Show top 20 intersections
      order.by = "freq",
      text.scale = 1.5,
      mainbar.y.label = "Proteins in Intersection",
      sets.x.label = "Total Proteins Annotated by Source")
dev.off() # Close PDF device


# === Plot 3: Top N Terms Bar Plots ===
cat("Generating... 3. Top N Bar Plots\n")

# Helper function: to clean, count, and plot
plot_top_terms <- function(data, column_name, separator, n = 15, title) {
  
  item_col <- sym(column_name) # Convert string to column name
  
  # Clean and count
  counts_df <- data %>%
    # --- 修复：明确使用 dplyr::select 来避免包冲突 ---
    dplyr::select(Protein_ID, {{ item_col }}) %>%
    filter({{ item_col }} != "-") %>%
    # Remove descriptions in parentheses, e.g., (Terpene synthase...)
    mutate(Term = str_replace_all({{ item_col }}, "\\(.*?\\)", "")) %>%
    # Split multiple entries
    separate_rows(Term, sep = separator) %>%
    # Remove GO term source, e.g., (InterPro)
    mutate(Term = str_replace_all(Term, "\\(.*?\\)", "")) %>%
    # Extract GO ID
    mutate(Term = str_extract(Term, "GO:\\d+|IPR\\d+|PF\\d+|PTHR\\d+")) %>%
    filter(!is.na(Term) & Term != "") %>%
    count(Term, sort = TRUE, name = "Count") %>%
    slice_max(order_by = Count, n = n)
  
  # Plot
  p <- ggplot(counts_df, aes(x = fct_reorder(Term, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = title, x = "Term ID", y = "Protein Count") +
    theme_minimal()
  
  return(p)
}

# 3a. Top 15 Pfam
p_pfam <- plot_top_terms(data, "Pfam_Annotation", ";", 15, "Top 15 Pfam Domains")
# --- Modify: Remove ggsave ---
# ggsave(file.path(output_dir, "3a_top15_pfam.pdf"), p_pfam, width = 10, height = 7)

# 3b. Top 15 InterPro
p_interpro <- plot_top_terms(data, "InterPro_Entries", ";", 15, "Top 15 InterPro Entries")
# --- Modify: Remove ggsave ---
# ggsave(file.path(output_dir, "3b_top15_interpro.pdf"), p_interpro, width = 10, height = 7)

# 3c. Top 15 GO
p_go <- plot_top_terms(data, "GO_Terms", ";", 15, "Top 15 GO Terms")
# --- Modify: Remove ggsave ---
# ggsave(file.path(output_dir, "3c_top15_go.pdf"), p_go, width = 10, height = 7)


# === Plot 4: GO Ontology Bar Plot ===
cat("Generating... 4. GO Ontology Bar Plot\n")
# First, get all unique GO terms and their total counts
all_go_counts <- data %>%
  dplyr::select(Protein_ID, GO_Terms) %>% # --- 修复：这里也明确使用 dplyr::select ---
  filter(GO_Terms != "-") %>%
  separate_rows(GO_Terms, sep = ";") %>%
  mutate(GO_ID = str_extract(GO_Terms, "GO:\\d+")) %>%
  filter(!is.na(GO_ID)) %>%
  count(GO_ID, sort = TRUE, name = "Count")

# Get list of GO IDs
go_ids <- all_go_counts$GO_ID
# Map GO IDs to their ontology (BP, CC, MF)
# keytype must be "GOID"
# --- 注意：这里我们 *特意* 使用 AnnotationDbi::select，这是正确的 ---
go_map <- AnnotationDbi::select(GO.db,
                                keys = go_ids,
                                columns = "ONTOLOGY",
                                keytype = "GOID")

# Join ontology info back to our counts
go_ontology_summary <- all_go_counts %>%
  left_join(go_map, by = c("GO_ID" = "GOID")) %>%
  filter(!is.na(ONTOLOGY)) %>% # Filter out potentially obsolete GO IDs
  # Summarize by major category
  group_by(ONTOLOGY) %>%
  # Note: We are summing the total occurrences of each term, not the count of unique terms
  summarise(TotalAnnotations = sum(Count)) %>%
  mutate(Category = case_when(
    ONTOLOGY == "BP" ~ "Biological Process (BP)",
    ONTOLOGY == "MF" ~ "Molecular Function (MF)",
    ONTOLOGY == "CC" ~ "Cellular Component (CC)",
    TRUE ~ ONTOLOGY
  ))

# Plot GO ontology breakdown
p4 <- ggplot(go_ontology_summary, aes(x = Category, y = TotalAnnotations, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "GO Annotation Ontology Breakdown",
       x = "GO Ontology",
       y = "Total Annotation Entries") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# --- Modify: Remove ggsave ---
# ggsave(file.path(output_dir, "4_go_ontology_breakdown.pdf"), p4, width = 8, height = 6)


# === NEW: Plot 5: Combine Plots (with patchwork) ===
cat("Generating... 5. Combining all ggplot plots\n")
# UpsetR plot (Plot 2) cannot be combined as it is not a ggplot object.
# We will combine p1, p4, p_pfam, p_interpro, p_go

# Define a layout:
# (p1) (p4)
# (p_pfam) (p_interpro)
# (p_go) (blank)
layout <- (p1 | p4) /
  (p_pfam | p_interpro) /
  (p_go | plot_spacer()) # plot_spacer() is a placeholder from patchwork

# Add a main title and auto-tagging (A, B, C...)
combined_plot <- layout + 
  plot_annotation(
    title = 'InterProScan Annotation Overview',
    tag_levels = 'A' # A, B, C, D, E...
  ) & 
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

# Save the combined plot
# Adjust size for a 3-row, 2-column layout
ggsave(file.path(output_dir, "5_combined_annotation_summary.pdf"), 
       combined_plot, 
       width = 16, 
       height = 18)
ggsave(file.path(output_dir, "5_combined_annotation_summary.png"), 
       combined_plot, 
       dpi = 1000,
       width = 16, 
       height = 18)

cat("All done!\n")