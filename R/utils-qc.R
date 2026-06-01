make_flow_tables <- function(vdj_list, seurat_filt, seurat_merged = NULL,
                             config_regex, prefix_regex, suffix_regex) {
  bcr_annotated <- vdj_list$bcr_annotated
  bcr_productive <- vdj_list$bcr_productive
  barcode_summary <- vdj_list$barcode_summary

  filtered_lookup <- make_gex_lookup(
    seurat_obj = seurat_filt,
    config_regex = config_regex,
    prefix_regex = prefix_regex,
    suffix_regex = suffix_regex,
    cell_col_name = "filtered_gex_cell"
  )

  merged_lookup <- make_gex_lookup(
    seurat_obj = seurat_merged,
    config_regex = config_regex,
    prefix_regex = prefix_regex,
    suffix_regex = suffix_regex,
    cell_col_name = "merged_gex_cell"
  )

  paired_vdj <- barcode_summary |>
    dplyr::filter(paired_bcr)

  vdj_gex_match <- paired_vdj |>
    dplyr::left_join(
      filtered_lookup,
      by = c("config", "barcode_clean")
    ) |>
    dplyr::left_join(
      merged_lookup,
      by = c("config", "barcode_clean")
    ) |>
    dplyr::mutate(
      maps_to_filtered_gex = !is.na(filtered_gex_cell),
      maps_to_merged_gex = !is.na(merged_gex_cell)
    )

  total_contigs <- nrow(bcr_annotated)
  productive_full_length_contigs <- nrow(bcr_productive)
  unique_vdj_barcodes <- nrow(barcode_summary)
  paired_vdj_barcodes <- sum(barcode_summary$paired_bcr)
  unpaired_vdj_barcodes <- unique_vdj_barcodes - paired_vdj_barcodes
  mapped_to_filtered_gex <- sum(vdj_gex_match$maps_to_filtered_gex)
  not_mapped_to_filtered_gex <- sum(!vdj_gex_match$maps_to_filtered_gex)
  total_filtered_gex_cells <- ncol(seurat_filt)
  total_merged_gex_cells <- ifelse(is.null(seurat_merged), NA_integer_, ncol(seurat_merged))
  mapped_to_merged_gex <- ifelse(is.null(seurat_merged), NA_integer_, sum(vdj_gex_match$maps_to_merged_gex))

  flow_summary <- tibble::tibble(
    metric = c(
      "Total VDJ contig rows",
      "Productive full-length contig rows",
      "Unique VDJ barcodes",
      "Paired heavy + light VDJ barcodes",
      "Unpaired/incomplete VDJ barcodes",
      "Total filtered GEX cells",
      "Total merged/pre-QC GEX cells",
      "Paired VDJ barcodes mapped to merged/pre-QC GEX",
      "Paired VDJ barcodes mapped to filtered GEX",
      "Paired VDJ barcodes not mapped to filtered GEX",
      "Percent paired VDJ mapped to filtered GEX"
    ),
    value = c(
      total_contigs,
      productive_full_length_contigs,
      unique_vdj_barcodes,
      paired_vdj_barcodes,
      unpaired_vdj_barcodes,
      total_filtered_gex_cells,
      total_merged_gex_cells,
      mapped_to_merged_gex,
      mapped_to_filtered_gex,
      not_mapped_to_filtered_gex,
      round(mapped_to_filtered_gex / paired_vdj_barcodes * 100, 2)
    )
  )

  barcode_loss_summary <- make_barcode_loss_summary(
    paired_vdj_barcodes = paired_vdj_barcodes,
    vdj_gex_match = vdj_gex_match,
    seurat_merged = seurat_merged
  )

  clone_table <- bcr_productive |>
    dplyr::filter(!is.na(raw_clonotype_id), raw_clonotype_id != "") |>
    dplyr::distinct(config, barcode, barcode_clean, raw_clonotype_id)

  clone_sizes <- clone_table |>
    dplyr::count(raw_clonotype_id, name = "clone_size") |>
    dplyr::arrange(dplyr::desc(clone_size))

  vdj_qc_summary <- tibble::tibble(
    metric = c(
      "Total contigs",
      "Unique barcodes",
      "Productive full-length contigs",
      "Unique productive full-length barcodes",
      "Paired heavy + light barcodes",
      "IGH contigs",
      "IGK contigs",
      "IGL contigs",
      "Barcodes with duplicate IGH",
      "Barcodes with duplicate light chain",
      "Barcodes with >2 productive contigs",
      "Unique clonotypes",
      "Expanded clonotypes",
      "Singleton clonotypes",
      "Largest clonotype size"
    ),
    value = c(
      nrow(bcr_annotated),
      dplyr::n_distinct(paste(bcr_annotated$config, bcr_annotated$barcode_clean, sep = "_")),
      nrow(bcr_productive),
      dplyr::n_distinct(paste(bcr_productive$config, bcr_productive$barcode_clean, sep = "_")),
      sum(barcode_summary$paired_bcr),
      sum(bcr_productive$chain == "IGH"),
      sum(bcr_productive$chain == "IGK"),
      sum(bcr_productive$chain == "IGL"),
      barcode_summary |> dplyr::filter(n_heavy > 1) |> nrow(),
      barcode_summary |> dplyr::filter(n_light > 1) |> nrow(),
      barcode_summary |> dplyr::filter(n_contigs > 2) |> nrow(),
      nrow(clone_sizes),
      clone_sizes |> dplyr::filter(clone_size > 1) |> nrow(),
      clone_sizes |> dplyr::filter(clone_size == 1) |> nrow(),
      ifelse(nrow(clone_sizes) > 0, max(clone_sizes$clone_size), NA_integer_)
    )
  )

  by_config <- barcode_summary |>
    dplyr::group_by(config) |>
    dplyr::summarise(
      unique_vdj_barcodes = dplyr::n(),
      paired_bcr_barcodes = sum(paired_bcr),
      unpaired_barcodes = sum(!paired_bcr),
      barcodes_with_duplicate_igh = sum(n_heavy > 1),
      barcodes_with_duplicate_light = sum(n_light > 1),
      barcodes_with_more_than_2_contigs = sum(n_contigs > 2),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      vdj_gex_match |>
        dplyr::group_by(config) |>
        dplyr::summarise(
          paired_vdj_mapped_to_filtered_gex = sum(maps_to_filtered_gex),
          paired_vdj_not_mapped_to_filtered_gex = sum(!maps_to_filtered_gex),
          percent_paired_mapped_to_filtered_gex = round(
            paired_vdj_mapped_to_filtered_gex / dplyr::n() * 100,
            2
          ),
          paired_vdj_mapped_to_merged_gex = sum(maps_to_merged_gex),
          .groups = "drop"
        ),
      by = "config"
    )

  v_gene_usage <- bcr_productive |>
    dplyr::filter(!is.na(v_gene), v_gene != "") |>
    dplyr::count(chain, v_gene, name = "n") |>
    dplyr::arrange(chain, dplyr::desc(n))

  j_gene_usage <- bcr_productive |>
    dplyr::filter(!is.na(j_gene), j_gene != "") |>
    dplyr::count(chain, j_gene, name = "n") |>
    dplyr::arrange(chain, dplyr::desc(n))

  list(
    flow_summary = flow_summary,
    barcode_loss_summary = barcode_loss_summary,
    vdj_qc_summary = vdj_qc_summary,
    by_config = by_config,
    vdj_gex_match = vdj_gex_match,
    filtered_lookup = filtered_lookup,
    merged_lookup = merged_lookup,
    barcode_summary = barcode_summary,
    clone_sizes = clone_sizes,
    v_gene_usage = v_gene_usage,
    j_gene_usage = j_gene_usage
  )
}

make_barcode_loss_summary <- function(paired_vdj_barcodes, vdj_gex_match, seurat_merged = NULL) {
  mapped_filtered <- sum(vdj_gex_match$maps_to_filtered_gex)
  not_mapped_filtered <- sum(!vdj_gex_match$maps_to_filtered_gex)

  if (is.null(seurat_merged)) {
    tibble::tibble(
      stage = c(
        "Paired VDJ barcodes",
        "Mapped to filtered/final GEX",
        "Not mapped to filtered/final GEX"
      ),
      count = c(
        paired_vdj_barcodes,
        mapped_filtered,
        not_mapped_filtered
      ),
      percent_of_paired_vdj = round(count / paired_vdj_barcodes * 100, 2)
    )
  } else {
    mapped_merged <- sum(vdj_gex_match$maps_to_merged_gex)
    in_merged_not_filtered <- sum(vdj_gex_match$maps_to_merged_gex & !vdj_gex_match$maps_to_filtered_gex)
    not_in_merged <- sum(!vdj_gex_match$maps_to_merged_gex)

    tibble::tibble(
      stage = c(
        "Paired VDJ barcodes",
        "Mapped to merged/pre-QC GEX",
        "Mapped to filtered/final GEX",
        "Present in merged but lost before filtered/final GEX",
        "Not found in merged/pre-QC GEX"
      ),
      count = c(
        paired_vdj_barcodes,
        mapped_merged,
        mapped_filtered,
        in_merged_not_filtered,
        not_in_merged
      ),
      percent_of_paired_vdj = round(count / paired_vdj_barcodes * 100, 2)
    )
  }
}
