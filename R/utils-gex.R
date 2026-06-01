make_gex_lookup <- function(seurat_obj, config_regex, prefix_regex, suffix_regex, cell_col_name = "gex_cell") {
  if (is.null(seurat_obj)) {
    return(tibble::tibble(
      !!cell_col_name := character(),
      config = character(),
      barcode_clean = character()
    ))
  }

  gex_cells <- colnames(seurat_obj)

  tibble::tibble(
    !!cell_col_name := gex_cells,
    config = stringr::str_extract(gex_cells, config_regex),
    barcode_clean = gex_cells |>
      stringr::str_remove(prefix_regex) |>
      stringr::str_remove(suffix_regex)
  )
}
