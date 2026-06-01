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
