#' Save BUSCO plots as PNG and/or PDF
#' @param p ggplot object.
#' @param out_prefix Output path prefix without extension.
#' @param width Figure width.
#' @param height Figure height.
#' @param dpi DPI for bitmap outputs.
#' @param save_png Logical.
#' @param save_pdf Logical.
#' @return Invisibly returns output file paths.
#' @export
save_busco_plots <- function(p, out_prefix, width = 7, height = 6.5, dpi = 1000, save_png = TRUE, save_pdf = TRUE) {
  pkgs <- c("ggplot2","fs")
  invisible(lapply(pkgs, requireNamespace, quietly = TRUE))
  fs::dir_create(fs::path_dir(out_prefix))

  out_files <- character()
  if (save_png) {
    pngf <- paste0(out_prefix, ".png")
    ggplot2::ggsave(pngf, plot = p, width = width, height = height, dpi = dpi, bg = "white")
    out_files <- c(out_files, pngf)
  }
  if (save_pdf) {
    pdff <- paste0(out_prefix, ".pdf")
    ggplot2::ggsave(pdff, plot = p, width = width, height = height, dpi = dpi)
    out_files <- c(out_files, pdff)
  }
  invisible(out_files)
}
