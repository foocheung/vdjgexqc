# =============================================================================
# utils-cite.R
# CITE-seq / ADT antigen scoring and HSA specificity calls
# Antigens are detected dynamically — no hardcoded DV1-DV4 names.
# =============================================================================

# ── Feature detection ─────────────────────────────────────────────────────────

#' Find ADT feature names matching a regex pattern
find_adt_features <- function(pattern, feature_names) {
  grep(pattern, feature_names, value = TRUE, ignore.case = TRUE, perl = TRUE)
}

#' Return all non-HSA ADT feature names from a Seurat object.
#' HSA is identified by common naming conventions; everything else is
#' treated as a potential antigen the user can select from.
detect_adt_features <- function(seurat_obj, assay = "ADT") {
  if (!assay %in% names(seurat_obj@assays)) return(NULL)
  feat <- rownames(seurat_obj[[assay]])
  
  hsa_features <- find_adt_features(
    "\\bHSA\\b|Human[._]Serum[._]Albumin|albumin", feat
  )
  antigen_features <- setdiff(feat, hsa_features)
  
  list(
    all_features     = feat,
    antigen_features = antigen_features,   # user selects from these
    hsa_features     = hsa_features
  )
}

# ── Per-cell scoring ──────────────────────────────────────────────────────────

#' Extract per-cell max signal for a set of ADT features.
#' Returns a numeric vector length == ncol(adt_mat).
get_antigen_score <- function(adt_mat, features) {
  features <- intersect(features, rownames(adt_mat))
  if (length(features) == 0) return(rep(NA_real_, ncol(adt_mat)))
  if (length(features) == 1) return(as.numeric(adt_mat[features, ]))
  apply(adt_mat[features, , drop = FALSE], 2, max, na.rm = TRUE)
}

#' Compute per-cell scores for each user-selected antigen plus HSA.
#' `selected_antigens`: character vector of individual ADT feature names chosen
#'   by the user (one score column per feature).
#' `hsa_features`: character vector of HSA feature name(s) — collapsed to max.
#' Returns a tibble: cell, one column per antigen, HSA, plus any meta_cols.
compute_antigen_scores <- function(seurat_obj,
                                   selected_antigens,
                                   hsa_features,
                                   assay     = "ADT",
                                   meta_cols = c("subject", "condition",
                                                 "monaco_main", "monaco_fine")) {
  adt_data <- tryCatch(
    SeuratObject::GetAssayData(seurat_obj, assay = assay, slot  = "data"),
    error = function(e)
      SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = "data")
  )
  adt_mat <- as.matrix(adt_data)
  
  # One column per selected antigen
  antigen_scores <- purrr::map_dfc(
    stats::setNames(selected_antigens, selected_antigens),
    ~ get_antigen_score(adt_mat, .x)
  )
  
  scores <- tibble::tibble(cell = colnames(seurat_obj)) |>
    dplyr::bind_cols(antigen_scores) |>
    dplyr::mutate(HSA = get_antigen_score(adt_mat, hsa_features))
  
  present_cols <- intersect(meta_cols, colnames(seurat_obj@meta.data))
  if (length(present_cols) > 0) {
    meta <- seurat_obj@meta.data[, present_cols, drop = FALSE] |>
      tibble::rownames_to_column("cell")
    scores <- dplyr::left_join(scores, meta, by = "cell")
  }
  scores
}

# ── Threshold grid and calling ────────────────────────────────────────────────

default_threshold_grid <- function() {
  tibble::tribble(
    ~threshold_label,       ~min_ag_minus_hsa, ~min_ag_to_hsa_ratio,
    "loose_0.5diff_1.5x",            0.5,               1.5,
    "medium_1diff_2x",               1.0,               2.0,
    "strict_1.5diff_3x",             1.5,               3.0,
    "very_strict_2diff_4x",          2.0,               4.0
  )
}

#' Pivot scores to long form and apply all thresholds.
#' `antigen_cols`: names of antigen columns in ag_scores (everything except
#'   cell, HSA, and metadata columns).
apply_antigen_thresholds <- function(ag_scores,
                                     antigen_cols,
                                     threshold_grid = default_threshold_grid()) {
  long <- ag_scores |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(antigen_cols),
      names_to  = "antigen",
      values_to = "ag_value"
    ) |>
    dplyr::mutate(
      ag_minus_HSA    = ag_value - HSA,
      ag_to_HSA_ratio = exp(ag_value) / pmax(exp(HSA), 1e-6)
    )
  
  tidyr::crossing(long, threshold_grid) |>
    dplyr::mutate(
      ag_specific = !is.na(ag_value) & !is.na(HSA) &
        ag_minus_HSA    >= min_ag_minus_hsa &
        ag_to_HSA_ratio >= min_ag_to_hsa_ratio
    )
}

#' Collapse per-antigen calls to per-cell "any antigen specific" calls.
collapse_antigen_calls <- function(ag_long) {
  ag_long |>
    dplyr::group_by(threshold_label, cell) |>
    dplyr::summarise(
      ag_any_specific      = any(ag_specific, na.rm = TRUE),
      ag_specific_antigens = paste(antigen[ag_specific], collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::mutate(ag_specific_antigens = dplyr::na_if(ag_specific_antigens, ""))
}

#' Write calls at a single threshold back onto Seurat metadata.
apply_antigen_calls_to_seurat <- function(seurat_obj, ag_calls_wide,
                                          chosen_threshold = "medium_1diff_2x") {
  calls <- ag_calls_wide |> dplyr::filter(threshold_label == chosen_threshold)
  idx   <- match(colnames(seurat_obj), calls$cell)
  
  seurat_obj$ag_threshold_used    <- chosen_threshold
  seurat_obj$ag_any_specific      <- calls$ag_any_specific[idx]
  seurat_obj$ag_specific_antigens <- calls$ag_specific_antigens[idx]
  seurat_obj$ag_positive          <- dplyr::if_else(
    seurat_obj$ag_any_specific == TRUE, "ag_pos", "ag_neg"
  )
  seurat_obj
}

# ── Summary tables ────────────────────────────────────────────────────────────

make_threshold_sensitivity <- function(ag_calls_wide, condition_col = "condition") {
  if (!condition_col %in% colnames(ag_calls_wide)) {
    ag_calls_wide[[condition_col]] <- "All cells"
  }
  ag_calls_wide |>
    dplyr::count(threshold_label,
                 condition = .data[[condition_col]],
                 ag_any_specific, name = "n") |>
    dplyr::group_by(threshold_label, condition) |>
    dplyr::mutate(pct = round(100 * n / sum(n), 2)) |>
    dplyr::ungroup() |>
    dplyr::filter(ag_any_specific)
}

make_ag_subject_summary <- function(seurat_obj,
                                    condition_filter = "Antigen positive") {
  meta <- seurat_obj@meta.data
  if ("condition" %in% colnames(meta))
    meta <- dplyr::filter(meta, condition == condition_filter)
  if (!"ag_positive" %in% colnames(meta)) return(NULL)
  meta |>
    tibble::rownames_to_column("cell") |>
    dplyr::filter(!is.na(subject), !is.na(ag_positive)) |>
    dplyr::count(subject, ag_positive, name = "n") |>
    dplyr::group_by(subject) |>
    dplyr::mutate(total_cells = sum(n), percent = round(100 * n / total_cells, 2)) |>
    dplyr::ungroup()
}

make_cross_reactivity_summary <- function(seurat_obj,
                                          condition_filter = "Antigen positive") {
  meta <- seurat_obj@meta.data |> tibble::rownames_to_column("cell")
  if ("condition" %in% colnames(meta))
    meta <- dplyr::filter(meta, condition == condition_filter)
  if (!"ag_specific_antigens" %in% colnames(meta)) return(NULL)
  meta |>
    dplyr::filter(ag_positive == "ag_pos", !is.na(ag_specific_antigens)) |>
    dplyr::mutate(
      n_antigens = stringr::str_count(ag_specific_antigens, ";") + 1,
      reactivity = dplyr::case_when(
        n_antigens == 1 ~ "Mono-reactive",
        n_antigens == 2 ~ "Bi-reactive",
        n_antigens >= 3 ~ "Multi-reactive (\u22653)"
      )
    ) |>
    dplyr::count(subject, reactivity, name = "n") |>
    dplyr::group_by(subject) |>
    dplyr::mutate(pct = round(100 * n / sum(n), 2)) |>
    dplyr::ungroup()
}

# ── B-cell subtype comparison ─────────────────────────────────────────────────

assign_comparison_groups <- function(seurat_obj,
                                     condition_col    = "condition",
                                     ag_pos_value     = NULL,   # user-supplied level
                                     ag_neg_value     = NULL) { # user-supplied level
  meta    <- seurat_obj@meta.data
  missing <- setdiff(c(condition_col, "ag_positive"), colnames(meta))
  if (length(missing) > 0) {
    warning("Cannot assign comparison groups — missing: ", paste(missing, collapse = ", "))
    return(seurat_obj)
  }
  
  cond <- meta[[condition_col]]
  
  seurat_obj$comparison_group <- dplyr::case_when(
    !is.null(ag_neg_value) & cond == ag_neg_value                              ~ "Ag_neg_sort",
    !is.null(ag_pos_value) & cond == ag_pos_value & meta$ag_positive == "ag_pos" ~ "Ag_pos",
    !is.null(ag_pos_value) & cond == ag_pos_value                              ~ "Ag_pos_ag_neg",
    TRUE ~ NA_character_
  )
  seurat_obj
}

make_bcell_freq <- function(seurat_obj,
                            groups      = c("Ag_pos", "Ag_neg_sort"),
                            bcell_types = c("Naive", "Non-switched Mem",
                                            "Switched Mem", "Exhausted")) {
  meta    <- seurat_obj@meta.data |> tibble::rownames_to_column("cell")
  missing <- setdiff(c("subject", "comparison_group", "bcell_type"), colnames(meta))
  if (length(missing) > 0) {
    warning("Missing columns for B-cell freq: ", paste(missing, collapse = ", "))
    return(NULL)
  }
  meta |>
    dplyr::filter(comparison_group %in% groups, bcell_type %in% bcell_types) |>
    dplyr::count(subject, comparison_group, bcell_type) |>
    dplyr::group_by(subject, comparison_group) |>
    dplyr::mutate(percent = round(100 * n / sum(n), 2)) |>
    dplyr::ungroup()
}

# ── Main pipeline entry point ─────────────────────────────────────────────────

#' Run the full CITE antigen scoring pipeline.
#' @param selected_antigens Character vector of ADT feature names chosen by the
#'   user (from detect_adt_features()$antigen_features).
run_cite_pipeline <- function(seurat_obj,
                              selected_antigens,
                              hsa_features_override = NULL,
                              chosen_threshold = "medium_1diff_2x",
                              adt_assay        = "ADT",
                              condition_col    = "condition",
                              ag_pos_value     = NULL,
                              ag_neg_value     = NULL) {
  adt_info <- detect_adt_features(seurat_obj, assay = adt_assay)
  if (is.null(adt_info))
    stop("ADT assay '", adt_assay, "' not found in Seurat object.")
  
  # Prefer user-specified control antigens; fall back to auto-detected HSA
  hsa_feats <- if (!is.null(hsa_features_override) && length(hsa_features_override) > 0)
    hsa_features_override else adt_info$hsa_features
  
  if (length(hsa_feats) == 0)
    stop("No control antigen (HSA) detected or selected.")
  if (length(selected_antigens) == 0)
    stop("No antigens selected.")
  
  meta_cols_want <- c("subject", "condition", "monaco_main", "monaco_fine")
  
  ag_scores <- compute_antigen_scores(
    seurat_obj        = seurat_obj,
    selected_antigens = selected_antigens,
    hsa_features      = hsa_feats,
    assay             = adt_assay,
    meta_cols         = meta_cols_want
  )
  
  threshold_grid <- default_threshold_grid()
  ag_long        <- apply_antigen_thresholds(ag_scores, selected_antigens, threshold_grid)
  ag_calls_wide  <- collapse_antigen_calls(ag_long)
  
  meta_slim     <- ag_scores |> dplyr::select(cell, dplyr::any_of(meta_cols_want))
  ag_calls_wide <- dplyr::left_join(ag_calls_wide, meta_slim, by = "cell")
  
  seurat_obj <- apply_antigen_calls_to_seurat(seurat_obj, ag_calls_wide, chosen_threshold)
  seurat_obj <- assign_comparison_groups(seurat_obj,
                                         condition_col = condition_col,
                                         ag_pos_value  = ag_pos_value,
                                         ag_neg_value  = ag_neg_value)
  
  bcell_monaco_map <- c(
    "Exhausted B cells"           = "Exhausted",
    "Naive B cells"               = "Naive",
    "Non-switched memory B cells" = "Non-switched Mem",
    "Switched memory B cells"     = "Switched Mem"
  )
  if ("monaco_fine" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$bcell_type <- factor(
      dplyr::recode(seurat_obj$monaco_fine, !!!bcell_monaco_map),
      levels = c("Naive", "Non-switched Mem", "Switched Mem", "Exhausted")
    )
  }
  
  threshold_sensitivity <- make_threshold_sensitivity(
    ag_calls_wide,
    condition_col = if ("condition" %in% colnames(ag_calls_wide)) "condition" else "All cells"
  )
  
  list(
    seurat_obj            = seurat_obj,
    adt_info              = adt_info,
    selected_antigens     = selected_antigens,
    ag_scores             = ag_scores,
    ag_long               = ag_long,
    ag_calls_wide         = ag_calls_wide,
    threshold_grid        = threshold_grid,
    threshold_sensitivity = threshold_sensitivity,
    ag_subject_summary    = make_ag_subject_summary(seurat_obj),
    cross_react           = make_cross_reactivity_summary(seurat_obj),
    bcell_freq            = tryCatch(make_bcell_freq(seurat_obj), error = function(e) NULL),
    chosen_threshold      = chosen_threshold
  )
}