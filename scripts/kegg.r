# 确保加载所有需要的包
perform_kegg_enrichment <- function(
    deg_data,           # 必需：原始的差异表达结果 (e.g., CK_L_deg)
    comparison_name,    # 必需：比较组名称 (e.g., "CK_L_10_vs_CK_L_3")
    direction,          # 必需："UP" 或 "DOWN"
    output_dir,         # 必需：保存文件的目录 (e.g., save_dir)
    term2gene_map,      # 必需：你的 term2gene 映射文件
    term2name_map,      # 必需：你的 term2name 映射文件
    gene2ko_map,        # 必需：你的 gene2ko 映射文件 (用于背景)
    
    padj_cutoff = 0.05,        # 可选：padj 阈值
    lfc_threshold = 1.0,       # 可选：log2FoldChange 阈值
    show_n_categories = 10,    # 可选：图表上展示的通路数量
    
    plot_width = 12,           # 可选：图片宽度
    plot_height = 6,           # 可选：图片高度
    plot_dpi = 1000            # 可选：图片 DPI
) {
  # loading packages
  require(dplyr)
  require(tibble)
  require(rlang)
  require(clusterProfiler)
  require(patchwork)

  direction_upper <- toupper(direction) 
  
  if (direction_upper == "UP") {
    # 使用 rlang::expr 来创建表达式，以便在 filter 中正确评估
    lfc_filter_expr <- rlang::expr(log2FoldChange > !!lfc_threshold)
    file_label <- "UP_GENE_KEGG"
    plot_title_label <- "UP GENE KEGG"
  } else if (direction_upper == "DOWN") {
    lfc_filter_expr <- rlang::expr(log2FoldChange < -!!lfc_threshold)
    file_label <- "DOWN_GENE_KEGG"
    plot_title_label <- "DOWN GENE KEGG"
  } else {
    stop("错误: 'direction' 参数必须是 'UP' 或 'DOWN'.")
  }
  
  # --- 2. 筛选差异基因 ---
  # (这里保留了你代码中特有的 GeneID -> KEGG_ID 转换逻辑)
  degs <- deg_data |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "GeneID") |>
    filter(padj < !!padj_cutoff & !!lfc_filter_expr) |>
    mutate(KEGG_ID = gsub("^.*_rb", "rb", GeneID)) |>
    pull(KEGG_ID)
  
  # 如果没有差异基因，则停止并给出提示
  if (length(degs) == 0) {
    print(paste("在", comparison_name, direction_upper, "中没有找到符合条件的基因。跳过分析。"))
    return(NULL)
  }
  
  # --- 3. 设置背景基因 ---
  background_genes <- unique(gene2ko_map$TranscriptID)
  
  # --- 4. 运行富集分析 ---
  print(paste("正在对", comparison_name, direction_upper, "的", length(degs), "个基因进行富集分析..."))
  enrich_result <- enricher(
    gene = degs,
    TERM2GENE = term2gene_map,
    TERM2NAME = term2name_map,
    pvalueCutoff = 0.05,  # (根据你的代码，这些是固定的)
    qvalueCutoff = 0.2,
    universe = background_genes
  )
  
  enrich_df <- as.data.frame(enrich_result)
  
  # --- 5. 保存CSV结果 ---
  file_base <- paste(comparison_name, file_label, sep = "_")
  csv_file <- file.path(output_dir, paste0(file_base, ".csv"))
  
  write.csv(enrich_df, row.names = FALSE, file = csv_file)
  print(paste("CSV 结果已保存至:", csv_file))
  
  # --- 6. 绘图和保存 ---
  # 检查是否有富集结果
  if (is.null(enrich_result) || nrow(enrich_df) == 0) {
    print("没有显著富集的结果，跳过绘图。")
    return(enrich_result) # 仍然返回空的结果对象
  }
  
  print("正在绘制图表...")
  p1 <- dotplot(enrich_result, showCategory = show_n_categories)
  p2 <- barplot(enrich_result, showCategory = show_n_categories)
  
  plot_title <- paste(comparison_name, plot_title_label, sep = " ")
  
  all_plots <- p1 + p2 +
    plot_layout() +
    plot_annotation(tag_levels = 'a', title = plot_title) &
    theme(plot.title = element_text(hjust = 0.5))
  
  # 保存 PNG
  png_file <- file.path(output_dir, paste0(file_base, ".png"))
  ggsave(png_file,
         width = plot_width, height = plot_height, dpi = plot_dpi, plot = all_plots)
  
  # 保存 PDF
  pdf_file <- file.path(output_dir, paste0(file_base, ".pdf"))
  ggsave(pdf_file,
         width = plot_width, height = plot_height, plot = all_plots)
  
  print(paste("图表 (PNG/PDF) 已保存至:", output_dir))
  
  # --- 7. 返回结果对象 ---
  # 返回富集分析的结果，以便在R环境中查看
  return(enrich_result)
}