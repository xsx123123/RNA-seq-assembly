#' Plot BUSCO composition across samples as alluvial diagram
#' @param df Wide data.frame from read_busco_jsons()
#' @param sample_order Optional character vector specifying sample display order.
#' @param category_order Order of BUSCO categories.
#' @param palette Fill colors for categories.
#' @param y_limit Numeric vector length 2, e.g., c(0,100).
#' @param label_percent Logical; show percentage labels on strata.
#' @return A ggplot object.
#' @export
plot_busco_alluvial <- function(
  df,
  sample_order = NULL,
  category_order = c("Missing", "Fragmented", "Multi.copy", "Single.copy"),
  palette = c("#ffb703","#fb8500","#8ecae6","#219ebc"),
  y_limit = c(0,100),
  label_percent = TRUE
) {
  pkgs <- c("dplyr","tidyr","stringr","ggplot2","ggalluvial","ggpubr","tibble","purrr")
  invisible(lapply(pkgs, requireNamespace, quietly = TRUE))

  # 长表
  long <- df |>
    dplyr::filter(rowname %in% c("Single.copy.percentage",
                                 "Multi.copy.percentage",
                                 "Fragmented.percentage",
                                 "Missing.percentage")) |>
    tidyr::pivot_longer(cols = -rowname, names_to = "Sample", values_to = "Percentage") |>
    dplyr::mutate(Category = stringr::str_remove(rowname, "\\.percentage$")) |>
    dplyr::mutate(Category = factor(Category, levels = category_order)) |>
    dplyr::mutate(Percentage = suppressWarnings(as.numeric(Percentage)))

  # 样本顺序
  if (!is.null(sample_order)) {
    long <- long |> dplyr::mutate(Sample = factor(Sample, levels = sample_order))
  } else {
    long <- long |> dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)))
  }

  p <- ggplot2::ggplot(
    data = long,
    ggplot2::aes(x = Sample, y = Percentage, fill = Category, alluvium = Category, stratum = Category)
  ) +
    ggalluvial::geom_alluvium(alpha = 0.7, color = "white", width = 0.4) +
    ggalluvial::geom_stratum(color = NA, width = 0.4) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(limits = y_limit, expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = "Percentage (%)", fill = "BUSCO Type") +
    ggpubr::theme_pubclean() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4))

  if (isTRUE(label_percent)) {
    p <- p + ggplot2::geom_text(
      stat = "stratum",
      ggplot2::aes(label = paste0(Percentage, "%")),
      check_overlap = TRUE, size = 3.5, color = "black"
    )
  }
  p
}
