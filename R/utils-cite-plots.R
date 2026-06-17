# =============================================================================
# utils-cite-plots.R
# Plotting functions reproducing the 303-6 CITE-seq slide figures
# =============================================================================

# ── Slide 1 / Image 1: Subject cell count table ───────────────────────────────

#' Table: n_cells_after_QC per subject (+ optional timepoint column)
make_subject_cell_table <- function(seurat_obj, subject_col, timepoint_col = NULL) {
  meta <- seurat_obj@meta.data |> tibble::rownames_to_column("cell")
  
  if (!subject_col %in% colnames(meta))
    stop("Subject column '", subject_col, "' not found in metadata.")
  
  grp_cols <- c(subject_col, timepoint_col)
  grp_cols <- grp_cols[!is.na(grp_cols) & grp_cols %in% colnames(meta)]
  
  meta |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) |>
    dplyr::summarise(n_cells_after_QC = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::across(dplyr::all_of(grp_cols[1])))
}

# ── Slide 3 / Image 3: ADT dot plot ──────────────────────────────────────────

#' Dot plot: average ADT expression × percent expressed per Monaco cell type
#' Reproduces Q1d from the slides.
make_adt_dotplot <- function(seurat_obj,
                             antigen_features,
                             celltype_col = "monaco_main",
                             adt_assay    = "ADT") {
  if (!adt_assay %in% names(seurat_obj@assays))
    stop("ADT assay '", adt_assay, "' not found.")
  if (!celltype_col %in% colnames(seurat_obj@meta.data))
    stop("Cell-type column '", celltype_col, "' not found in metadata.")
  
  adt_data <- tryCatch(
    SeuratObject::GetAssayData(seurat_obj, assay = adt_assay, slot  = "data"),
    error = function(e)
      SeuratObject::GetAssayData(seurat_obj, assay = adt_assay, layer = "data")
  )
  
  feats <- intersect(antigen_features, rownames(adt_data))
  if (length(feats) == 0) stop("None of the selected antigen features found in ADT assay.")
  
  mat   <- as.matrix(adt_data[feats, , drop = FALSE])
  meta  <- seurat_obj@meta.data |> tibble::rownames_to_column("cell")
  
  long <- mat |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("cell") |>
    tidyr::pivot_longer(-cell, names_to = "feature", values_to = "expr") |>
    dplyr::left_join(meta |> dplyr::select(cell, celltype = dplyr::all_of(celltype_col)),
                     by = "cell") |>
    dplyr::filter(!is.na(celltype))
  
  dot_df <- long |>
    dplyr::group_by(feature, celltype) |>
    dplyr::summarise(
      avg_expr     = mean(expr, na.rm = TRUE),
      pct_expressed = 100 * mean(expr > 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(feature) |>
    dplyr::mutate(avg_expr_scaled = scale(avg_expr)[, 1]) |>
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(
    dot_df,
    ggplot2::aes(x = celltype, y = feature,
                 size = pct_expressed, colour = avg_expr_scaled)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient2(low = "white", mid = "purple", high = "darkblue",
                                    midpoint = 0, name = "Average\nExpression") +
    ggplot2::scale_size_continuous(range = c(1, 8), name = "Percent\nExpressed") +
    ggplot2::labs(x = "Monaco cell type", y = "ADT surface marker",
                  title = "ADT validation of GEX-defined PBMC cell populations") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  plotly::ggplotly(p)
}

# ── Slide 4 / Image 4: B-cell subpopulation stacked bar ──────────────────────

#' Stacked bar: B-cell subtype % of all PBMCs per subject
make_bcell_abundance_plot <- function(seurat_obj,
                                      subject_col   = "subject",
                                      timepoint_col = NULL,
                                      bcell_col     = "monaco_fine",
                                      bcell_types   = c("Naive B cells",
                                                        "Non-switched memory B cells",
                                                        "Switched memory B cells",
                                                        "Exhausted B cells",
                                                        "Plasmablasts")) {
  meta      <- seurat_obj@meta.data |> tibble::rownames_to_column("cell")
  total_per <- meta |>
    dplyr::count(.data[[subject_col]], name = "total_cells")
  
  x_label <- if (!is.null(timepoint_col) && timepoint_col %in% colnames(meta)) {
    paste0(subject_col, "\n", timepoint_col)
  } else subject_col
  
  plot_df <- meta |>
    dplyr::filter(.data[[bcell_col]] %in% bcell_types) |>
    dplyr::count(.data[[subject_col]], .data[[bcell_col]], name = "n") |>
    dplyr::left_join(total_per, by = subject_col) |>
    dplyr::mutate(pct = 100 * n / total_cells)
  
  if (!is.null(timepoint_col) && timepoint_col %in% colnames(meta)) {
    tp_map <- meta |>
      dplyr::distinct(.data[[subject_col]], .data[[timepoint_col]]) |>
      dplyr::mutate(x_label = paste0(.data[[subject_col]], "\n", .data[[timepoint_col]]))
    plot_df <- dplyr::left_join(plot_df, tp_map, by = subject_col)
  } else {
    plot_df$x_label <- as.character(plot_df[[subject_col]])
  }
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = x_label, y = pct, fill = .data[[bcell_col]])
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(x = NULL, y = "% of all PBMCs",
                  fill = "B-cell subtype",
                  title = "B cell subpopulation abundance per subject") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  plotly::ggplotly(p)
}

# ── Slide 5 / Image 5: Normalised cell-type proportion heatmap ───────────────

#' Heatmap: z-scored cell-type frequencies per subject (clustered)
make_celltype_heatmap <- function(seurat_obj,
                                  subject_col  = "subject",
                                  celltype_col = "monaco_fine") {
  meta <- seurat_obj@meta.data |> tibble::rownames_to_column("cell")
  
  freq <- meta |>
    dplyr::filter(!is.na(.data[[celltype_col]]), !is.na(.data[[subject_col]])) |>
    dplyr::count(.data[[subject_col]], .data[[celltype_col]], name = "n") |>
    dplyr::group_by(.data[[subject_col]]) |>
    dplyr::mutate(pct = 100 * n / sum(n)) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(id_cols = subject_col,
                       names_from  = celltype_col,
                       values_from = "pct",
                       values_fill = 0)
  
  mat <- freq |>
    tibble::column_to_rownames(subject_col) |>
    as.matrix()
  
  # z-score per cell type (column)
  mat_z <- scale(mat)
  
  heatmaply::heatmaply(
    mat_z,
    colors       = colorRampPalette(c("blue", "white", "red"))(256),
    xlab         = "",
    ylab         = "",
    main         = "Normalised cell-type proportions per subject",
    show_dendrogram = c(TRUE, TRUE),
    key.title    = "Z-score"
  )
}

# ── Slide 7 / Image 7: Antigen vs HSA scatter grid ───────────────────────────

#' Facet scatter: each antigen signal vs HSA per subject
#' Red = above threshold (ag_specific == TRUE); grey = below.
make_ag_vs_hsa_plot <- function(ag_scores,
                                antigen_cols,
                                subject_col  = "subject",
                                ag_long      = NULL,
                                chosen_threshold = "medium_1diff_2x") {
  # Colour by specificity if ag_long provided
  if (!is.null(ag_long)) {
    spec <- ag_long |>
      dplyr::filter(threshold_label == chosen_threshold) |>
      dplyr::select(cell, antigen, ag_specific)
  }
  
  long <- ag_scores |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(antigen_cols),
      names_to  = "antigen",
      values_to = "ag_value"
    ) |>
    dplyr::filter(!is.na(ag_value), !is.na(HSA))
  
  if (!is.null(ag_long)) {
    long <- dplyr::left_join(long, spec, by = c("cell", "antigen"))
  } else {
    long$ag_specific <- FALSE
  }
  
  subj_col_sym <- rlang::sym(subject_col)
  
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = HSA, y = ag_value,
                 colour = ag_specific,
                 text   = paste0("Cell: ", cell,
                                 "<br>HSA: ", round(HSA, 2),
                                 "<br>", antigen, ": ", round(ag_value, 2)))
  ) +
    ggplot2::geom_point(size = 0.6, alpha = 0.6) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
    ggplot2::scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#CC0000"),
                                 guide  = "none") +
    ggplot2::facet_grid(
      rows = ggplot2::vars(antigen),
      cols = ggplot2::vars(!!subj_col_sym)
    ) +
    ggplot2::labs(x = "HSA signal, negative-control antigen",
                  y = "Antigen signal",
                  title = "Per-cell antigen vs HSA signal by subject") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
  
  plotly::ggplotly(p, tooltip = "text")
}

# ── Slide 8 / Image 8: B-cell subtype by comparison group ────────────────────

#' Stacked bar: B-cell subtype % per subject × comparison group (Ag_neg / DV_pos)
make_bcell_by_group_plot <- function(seurat_obj,
                                     subject_col = "subject",
                                     bcell_col   = "bcell_type",
                                     groups      = c("Ag_neg_sort", "Ag_pos"),
                                     bcell_types = c("Naive", "Non-switched Mem",
                                                     "Switched Mem", "Exhausted")) {
  meta <- seurat_obj@meta.data |> tibble::rownames_to_column("cell")
  
  required <- c(subject_col, "comparison_group", bcell_col)
  missing  <- setdiff(required, colnames(meta))
  if (length(missing) > 0)
    stop("Missing metadata columns: ", paste(missing, collapse = ", "))
  
  plot_df <- meta |>
    dplyr::filter(
      comparison_group %in% groups,
      .data[[bcell_col]] %in% bcell_types
    ) |>
    dplyr::count(.data[[subject_col]], comparison_group, .data[[bcell_col]], name = "n") |>
    dplyr::group_by(.data[[subject_col]], comparison_group) |>
    dplyr::mutate(pct = 100 * n / sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      comparison_group = factor(comparison_group, levels = groups)
    )
  
  colours <- c(
    "Naive"            = "#4472C4",
    "Non-switched Mem" = "#ED7D31",
    "Switched Mem"     = "#C00000",
    "Exhausted"        = "#70AD47"
  )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = comparison_group, y = pct,
                 fill = .data[[bcell_col]])
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = colours, name = "B-cell subtype") +
    ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(subject_col)), nrow = 2) +
    ggplot2::labs(x = NULL, y = "% of B cells",
                  title = "B-cell subtype by antigen-specificity group per subject") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  
  plotly::ggplotly(p)
}

# ── Slide 9 / Image 9: Top IGH clonotype lollipop ────────────────────────────

#' Lollipop: top IGH clonotypes in Ag-neg vs Ag-pos lanes for one subject.
#'
#' @param bcr_productive  Productive VDJ tibble from build_vdj_summary().
#' @param seurat_meta     Seurat metadata data.frame (rownames = cell barcodes),
#'                        must contain subject_col and comparison_group columns.
#' @param subject_id      The subject value to filter to.
#' @param subject_col     Name of the subject column in seurat_meta.
#' @param ag_pos_label    comparison_group value for antigen-positive cells.
#' @param ag_neg_label    comparison_group value for antigen-negative cells.
#' @param top_n           How many top clonotypes to show per lane.
make_clonotype_lollipop <- function(bcr_productive,
                                    seurat_meta,
                                    subject_id,
                                    subject_col   = "subject",
                                    ag_pos_label  = "Ag_pos",
                                    ag_neg_label  = "Ag_neg_sort",
                                    top_n         = 20) {
  
  # ── Resolve CDR3 column name ───────────────────────────────────────────────
  if (!"cdr3" %in% colnames(bcr_productive)) {
    cdr3_alt <- intersect(c("cdr3_aa", "junction_aa"), colnames(bcr_productive))
    if (length(cdr3_alt) == 0)
      stop("No CDR3 column found in VDJ data (expected 'cdr3', 'cdr3_aa', or 'junction_aa').")
    bcr_productive <- dplyr::rename(bcr_productive, cdr3 = dplyr::all_of(cdr3_alt[1]))
  }
  
  # ── Pull comparison_group + subject from Seurat metadata ──────────────────
  required_meta <- c(subject_col, "comparison_group")
  missing <- setdiff(required_meta, colnames(seurat_meta))
  if (length(missing) > 0)
    stop("Seurat metadata missing columns: ", paste(missing, collapse = ", "),
         ". Run the CITE pipeline first to assign comparison groups.")
  
  cell_groups <- seurat_meta |>
    tibble::rownames_to_column("cell") |>
    dplyr::filter(
      .data[[subject_col]] == subject_id,
      comparison_group %in% c(ag_pos_label, ag_neg_label)
    ) |>
    dplyr::select(cell, comparison_group) |>
    dplyr::mutate(
      # Extract the raw 16-base barcode (ACGT only), dropping any prefix/suffix
      join_key = stringr::str_extract(cell, "[ACGT]{16}")
    )
  
  if (nrow(cell_groups) == 0)
    stop("No cells found for subject '", subject_id,
         "' with comparison_group in c('", ag_pos_label, "', '", ag_neg_label, "').")
  
  # ── Subset VDJ to IGH, join on cleaned barcode ────────────────────────────
  igh_all <- bcr_productive |>
    dplyr::filter(chain == "IGH", !is.na(cdr3), !is.na(v_gene))
  
  igh <- igh_all |>
    dplyr::mutate(join_key = stringr::str_extract(barcode_clean, "[ACGT]{16}")) |>
    dplyr::inner_join(cell_groups, by = "join_key") |>
    dplyr::mutate(
      lane            = dplyr::if_else(comparison_group == ag_pos_label,
                                       "Ag-positive", "Ag-negative"),
      clonotype_label = paste0(cdr3, " | ", v_gene)
    )
  
  if (nrow(igh) == 0)
    stop(
      "No IGH contigs found for subject '", subject_id, "'.\n",
      "Sample VDJ barcode_clean: ", paste(head(igh_all$barcode_clean, 5), collapse = ", "), "\n",
      "Sample Seurat join_key:   ", paste(head(cell_groups$join_key,   5), collapse = ", "), "\n",
      "If these don't match, check the prefix/suffix regex settings."
    )
  
  # ── Top clonotypes independently per lane ─────────────────────────────────
  top_pos <- igh |> dplyr::filter(lane == "Ag-positive") |>
    dplyr::count(clonotype_label, name = "n") |>
    dplyr::slice_max(n, n = top_n, with_ties = FALSE) |>
    dplyr::pull(clonotype_label)
  
  top_neg <- igh |> dplyr::filter(lane == "Ag-negative") |>
    dplyr::count(clonotype_label, name = "n") |>
    dplyr::slice_max(n, n = top_n, with_ties = FALSE) |>
    dplyr::pull(clonotype_label)
  
  all_top <- union(top_pos, top_neg)
  
  # ── Build plot data ────────────────────────────────────────────────────────
  plot_df <- igh |>
    dplyr::filter(clonotype_label %in% all_top) |>
    dplyr::count(clonotype_label, lane, name = "n") |>
    tidyr::complete(clonotype_label, lane = c("Ag-positive", "Ag-negative"),
                    fill = list(n = 0)) |>
    dplyr::mutate(
      detected_in = dplyr::case_when(
        clonotype_label %in% top_pos & clonotype_label %in% top_neg ~ "Both lanes",
        clonotype_label %in% top_pos                                 ~ "Ag-pos only",
        TRUE                                                         ~ "Ag-neg only"
      ),
      # Order y-axis by total count across both lanes
      total = ave(n, clonotype_label, FUN = sum),
      clonotype_label = factor(clonotype_label,
                               levels = unique(clonotype_label[order(total)]))
    )
  
  colours <- c("Both lanes" = "black", "Ag-pos only" = "#00B0F0", "Ag-neg only" = "#FF6666")
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = lane, y = clonotype_label, size = n, colour = detected_in,
                 text = paste0(clonotype_label, "<br>Lane: ", lane, "<br>n = ", n))
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(range = c(2, 10), name = "Cell count") +
    ggplot2::scale_colour_manual(values = colours, name = "Detection status") +
    ggplot2::labs(x = NULL, y = "IGH CDR3 | V gene",
                  title = paste0("Top ", top_n, " IGH clonotypes per lane — subject: ", subject_id)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
  
  plotly::ggplotly(p, tooltip = "text")
}