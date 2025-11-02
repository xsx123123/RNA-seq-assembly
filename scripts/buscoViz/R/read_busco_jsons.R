#' Read multiple BUSCO JSON summaries and merge to wide format
#' @param files Character vector of JSON file paths.
#' @param sample_names Optional character vector of sample names; must match length(files).
#' @return A data.frame: columns are samples, rowname column stores BUSCO keys.
#' @export
read_busco_jsons <- function(files, sample_names = NULL) {
  stopifnot(length(files) >= 1)
  pkgs <- c("jsonlite","dplyr","tibble","purrr")
  invisible(lapply(pkgs, requireNamespace, quietly = TRUE))

  if (is.null(sample_names)) {
    sample_names <- basename(files)
  } else {
    stopifnot(length(sample_names) == length(files))
  }

  dfs <- purrr::map2(files, sample_names, function(f, nm) {
    js <- jsonlite::fromJSON(f)
    if (!"results" %in% names(js)) {
      cli::cli_abort("JSON does not contain 'results': {f}")
    }
    as.data.frame(js$results) |>
      t() |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "rowname") |>
      dplyr::rename(!!nm := V1)
  })

  df <- purrr::reduce(dfs, function(x, y) dplyr::left_join(x, y, by = "rowname"))
  df
}
