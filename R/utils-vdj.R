safe_bool <- function(x) {
  x %in% c(TRUE, "TRUE", "True", "true", "T", "t", "1", 1)
}

parse_config_from_file <- function(file_name) {
  cfg <- stringr::str_extract(basename(file_name), "(?<=multi_config_)[0-9]+")
  if (is.na(cfg)) cfg <- stringr::str_extract(basename(file_name), "(?<=config_)[0-9]+")
  if (is.na(cfg)) cfg <- stringr::str_extract(basename(file_name), "(?<=lane_)[0-9]+")
  if (is.na(cfg)) cfg <- stringr::str_extract(basename(file_name), "(?<=sample_)[0-9]+")
  cfg
}


# Derive which config number(s) correspond to the ag-pos sort lane.
#
# Strategy (in priority order):
#  1. If the Seurat condition column contains cells whose value matches
#     ag_pos_value, extract the config from those cell names — most precise.
#  2. Scan VDJ file names for keywords like "ag_pos", "agpos", "sorted",
#     "antigen_pos", etc. — handles descriptive file names.
#  3. Fall back to the highest numeric config in the uploaded file list —
#     matches the original convention (lane 7 = last uploaded file).
#
# Returns a character vector of config values (e.g. c("7")) that can be
# compared directly against the `config` column in bcr_annotated /
# barcode_summary / filtered_lookup.
derive_ag_pos_configs <- function(vdj_list, file_names, all_file_configs,
                                  ag_pos_value = NULL,
                                  seurat_filt  = NULL,
                                  condition_col = NULL) {
  
  bcr_annotated <- vdj_list$bcr_annotated
  
  # ── Strategy 1: match via Seurat condition column ────────────────────────
  if (!is.null(seurat_filt) && !is.null(condition_col) &&
      condition_col %in% colnames(seurat_filt@meta.data) &&
      !is.null(ag_pos_value) && nzchar(ag_pos_value)) {
    
    ag_pos_cells <- rownames(seurat_filt@meta.data)[
      !is.na(seurat_filt@meta.data[[condition_col]]) &
        seurat_filt@meta.data[[condition_col]] == ag_pos_value
    ]
    
    if (length(ag_pos_cells) > 0) {
      # Extract config numbers from these cell names using the same regex
      # parse_config_from_file uses, then intersect with known configs.
      cell_configs <- unique(stats::na.omit(
        stringr::str_extract(ag_pos_cells, "(?<=multi_config_)[0-9]+|(?<=config_)[0-9]+|(?<=lane_)[0-9]+|(?<=sample_)[0-9]+")
      ))
      valid <- intersect(cell_configs, unique(bcr_annotated$config))
      if (length(valid) > 0) {
        message("derive_ag_pos_configs: inferred ag-pos config(s) from Seurat condition column '",
                condition_col, "' == '", ag_pos_value, "': ",
                paste(valid, collapse = ", "))
        return(valid)
      }
    }
  }
  
  # ── Strategy 2: keyword scan of uploaded VDJ file names ──────────────────
  ag_pos_pattern <- stringr::regex(
    "ag.?pos|antigen.?pos|sorted|sort.?pos|agpos|pos.?sort",
    ignore_case = TRUE
  )
  keyword_hits <- stringr::str_detect(file_names, ag_pos_pattern)
  if (any(keyword_hits)) {
    configs_from_names <- unique(stats::na.omit(all_file_configs[keyword_hits]))
    valid <- intersect(configs_from_names, unique(bcr_annotated$config))
    if (length(valid) > 0) {
      message("derive_ag_pos_configs: inferred ag-pos config(s) from file name keywords: ",
              paste(valid, collapse = ", "),
              " (files: ", paste(file_names[keyword_hits], collapse = ", "), ")")
      return(valid)
    }
  }
  
  # ── Strategy 3: highest numeric config (original convention) ─────────────
  numeric_configs <- suppressWarnings(as.numeric(unique(stats::na.omit(all_file_configs))))
  if (length(numeric_configs) > 0) {
    highest <- as.character(max(numeric_configs))
    message("derive_ag_pos_configs: falling back to highest config number: ", highest)
    return(highest)
  }
  
  # Nothing worked — return NULL and let the Sankey functions handle it.
  message("derive_ag_pos_configs: could not determine ag-pos config; Sankey lane step will be skipped.")
  NULL
}

clean_vdj_barcode <- function(x) {
  x |>
    stringr::str_remove("-1$") |>
    stringr::str_remove("_1$")
}

build_vdj_summary <- function(vdj_files_df) {
  bcr_annotated <- purrr::map2_dfr(
    vdj_files_df$datapath,
    vdj_files_df$name,
    function(path, original_name) {
      readr::read_csv(path, show_col_types = FALSE) |>
        dplyr::mutate(
          source_file = original_name,
          config = parse_config_from_file(original_name)
        )
    }
  )
  
  required_cols <- c("barcode", "chain", "productive", "full_length")
  missing_cols <- setdiff(required_cols, colnames(bcr_annotated))
  
  if (length(missing_cols) > 0) {
    stop(
      paste(
        "VDJ file is missing required columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  file_config_lookup <- vdj_files_df |>
    dplyr::mutate(
      parsed_config = purrr::map_chr(name, parse_config_from_file),
      fallback_config = as.character(dplyr::row_number()),
      final_config = dplyr::if_else(is.na(parsed_config), fallback_config, parsed_config)
    ) |>
    dplyr::select(source_file = name, final_config)
  
  bcr_annotated <- bcr_annotated |>
    dplyr::left_join(file_config_lookup, by = "source_file") |>
    dplyr::mutate(
      config = dplyr::if_else(is.na(config), final_config, config),
      barcode_clean = clean_vdj_barcode(barcode)
    ) |>
    dplyr::select(-final_config)
  
  bcr_productive <- bcr_annotated |>
    dplyr::filter(
      safe_bool(productive),
      safe_bool(full_length)
    )
  
  for (missing_col in c("raw_clonotype_id", "v_gene", "j_gene")) {
    if (!missing_col %in% colnames(bcr_productive)) {
      bcr_productive[[missing_col]] <- NA_character_
    }
  }
  
  barcode_summary <- bcr_productive |>
    dplyr::group_by(config, barcode, barcode_clean) |>
    dplyr::summarise(
      has_heavy = any(chain == "IGH"),
      has_light = any(chain %in% c("IGK", "IGL")),
      n_contigs = dplyr::n(),
      n_heavy = sum(chain == "IGH"),
      n_light = sum(chain %in% c("IGK", "IGL")),
      clonotype_id = paste(unique(stats::na.omit(raw_clonotype_id)), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      paired_bcr = has_heavy & has_light,
      vdj_cell_id = paste0("multi_config_", config, "_", barcode_clean)
    )
  
  list(
    bcr_annotated = bcr_annotated,
    bcr_productive = bcr_productive,
    barcode_summary = barcode_summary
  )
}