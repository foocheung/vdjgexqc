make_gex_lookup <- function(seurat_obj, config_regex, prefix_regex, suffix_regex,
                            cell_col_name = "gex_cell") {
  if (is.null(seurat_obj)) {
    return(tibble::tibble(
      !!cell_col_name := character(),
      config           = character(),
      barcode_clean    = character()
    ))
  }
  
  gex_cells <- colnames(seurat_obj)
  
  # ── Extract config ──────────────────────────────────────────────────────────
  config_extracted <- stringr::str_extract(gex_cells, config_regex)
  
  # ── Extract barcode_clean ───────────────────────────────────────────────────
  # Strategy 1: apply user-supplied prefix/suffix regexes
  bc_regex <- gex_cells |>
    stringr::str_remove(prefix_regex) |>
    stringr::str_remove(suffix_regex)
  
  # Strategy 2: extract the raw 16-base [ACGT]{16} barcode as a robust fallback.
  # Used when the regex approach leaves non-barcode characters (e.g. lane prefix
  # like "AG_neg_" wasn't covered by the regex).
  bc_extract <- stringr::str_extract(gex_cells, "[ACGT]{16,18}")
  
  # Choose per-cell: use regex result if it looks like a pure barcode
  # (only A/C/G/T and optional trailing -digit), otherwise use extraction.
  is_clean <- stringr::str_detect(bc_regex, "^[ACGT]{14,18}(-[0-9]+)?$")
  barcode_clean <- dplyr::if_else(is_clean, bc_regex, bc_extract)
  
  tibble::tibble(
    !!cell_col_name := gex_cells,
    config        = config_extracted,
    barcode_clean = barcode_clean
  )
}