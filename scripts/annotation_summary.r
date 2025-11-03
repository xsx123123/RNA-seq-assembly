# --------------------------------------------------
# Step 1: Load required R packages
# --------------------------------------------------
# Required installations:
# install.packages("tidyverse")
# install.packages("patchwork")
# install.packages("stringr")
library(tidyverse)
library(patchwork)
library(stringr)

# --------------------------------------------------
# Step 2: Prepare and preprocess data
# --------------------------------------------------
root_dir <- '/data/jzhang/project/RNA-seq-assembly_2025.8.25/05.transcript_annotation'
data_path <- '/data/jzhang/project/RNA-seq-assembly_2025.8.25/05.transcript_annotation/TD2_pep_matches_annotated.tsv'

# Read and preprocess annotation data
data <- readr::read_tsv(data_path) %>%
  mutate(
    e_value = as.numeric(e_value),
    bit_score = as.numeric(bit_score),
    e_value_plot = if_else(e_value == 0, 1e-200, e_value)  # Handle zero e-values
  )

print("Data loaded and preprocessed.")

# --------------------------------------------------
# Plot 1: Top 10 Annotated Organisms
# --------------------------------------------------
top_organisms <- data %>%
  filter(!is.na(Organism)) %>%
  count(Organism, sort = TRUE) %>%
  top_n(10, n)

p_organisms <- ggplot(top_organisms, aes(x = reorder(Organism, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 Annotated Organisms",
    x = "Organism",
    y = "Gene Count"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

# --------------------------------------------------
# Plot 2: E-value Distribution
# --------------------------------------------------
p_evalue <- ggplot(data, aes(x = -log10(e_value_plot))) +
  geom_histogram(bins = 30, fill = "darkorange", alpha = 0.7, color = "white") +
  labs(
    title = "E-value Distribution",
    x = "-log10(E-value)",
    y = "Frequency"
  ) +
  geom_vline(xintercept = -log10(1e-10), linetype = "dashed", color = "red") +
  annotate("text", x = -log10(1e-10) + 5, y = Inf, 
           label = "1e-10", vjust = 2, color = "red", size = 3) +
  theme_minimal()

# --------------------------------------------------
# Plot 3: Bit Score Distribution
# --------------------------------------------------
p_bitscore <- ggplot(data, aes(x = bit_score)) +
  geom_histogram(bins = 30, fill = "forestgreen", alpha = 0.7, color = "white") +
  labs(
    title = "Bit Score Distribution",
    x = "Bit Score",
    y = "Frequency"
  ) +
  theme_minimal()

# --------------------------------------------------
# GO Annotation Visualization Function
# --------------------------------------------------
plot_go_terms <- function(data, col_name, plot_title, plot_color) {
  if (!col_name %in% names(data)) {
    warning(paste("Column", col_name, "not found. Skipping plot."))
    return(ggplot() + 
             labs(title = plot_title) + 
             annotate("text", x = 0.5, y = 0.5, label = "Data unavailable", size = 5) +
             theme_void())
  }
  
  go_counts <- data %>%
    select(!!sym(col_name)) %>%
    filter(!is.na(!!sym(col_name))) %>%
    separate_rows(!!sym(col_name), sep = "; ") %>%
    mutate(go_term = str_replace(!!sym(col_name), " \\[GO:\\d+\\]", "")) %>%
    filter(!is.na(go_term) & go_term != "NA" & go_term != "") %>%
    count(go_term, sort = TRUE) %>%
    top_n(10, n)
  
  if (nrow(go_counts) == 0) {
    return(ggplot() + 
             labs(title = plot_title) + 
             annotate("text", x = 0.5, y = 0.5, label = "No valid terms", size = 5) +
             theme_void())
  }
  
  ggplot(go_counts, aes(x = reorder(go_term, n), y = n)) +
    geom_bar(stat = "identity", fill = plot_color) +
    coord_flip() +
    labs(
      title = plot_title,
      x = "GO Term",
      y = "Count"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
}

# Generate GO plots
p_go_bp <- plot_go_terms(data, 
                         "Gene Ontology (biological process)", 
                         "Top 10 GO - Biological Process",
                         "purple")

p_go_cc <- plot_go_terms(data, 
                         "Gene Ontology (cellular component)", 
                         "Top 10 GO - Cellular Component",
                         "firebrick")

p_go_mf <- plot_go_terms(data, 
                         "Gene Ontology (molecular function)", 
                         "Top 10 GO - Molecular Function",
                         "darkcyan")

# --------------------------------------------------
# Step 3: Compose final visualization
# --------------------------------------------------
# Arrange plots in 3x2 grid layout
final_plot <- (p_organisms | p_evalue) /
               (p_go_bp | p_go_cc) /
               (p_go_mf | p_bitscore) +
  plot_annotation(
    title = "Transcriptome Annotation Summary",
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )

# Display and export results
print(final_plot)
ggsave(
  filename = file.path(root_dir, "annotation_summary.png"),
  plot = final_plot,
  width = 16,
  height = 12,
  dpi = 1000
)
ggsave(
  filename = file.path(root_dir, "annotation_summary.pdf"),
  plot = final_plot,
  width = 16,
  height = 12,
  dpi = 1000
)

print("Plot assembly complete. Results saved to annotation_summary.png and annotation_summary.pdf")
