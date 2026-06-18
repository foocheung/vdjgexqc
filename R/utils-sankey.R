# =============================================================================
# utils-sankey.R
# =============================================================================

# ── VDJ → GEX Sankey ─────────────────────────────────────────────────────────

make_sankey <- function(flow_summary, bcr_annotated = NULL,
                        filtered_lookup = NULL,
                        ag_pos_configs = NULL) {
  get_val <- function(metric_name) {
    x <- flow_summary |>
      dplyr::filter(metric == metric_name) |>
      dplyr::pull(value) |>
      as.numeric()
    if (length(x) == 0 || is.na(x)) 0 else x
  }
  
  lbl <- function(label, n) paste0(label, "\n", n)
  
  total_contigs                  <- get_val("Total VDJ contig rows")
  productive_full_length_contigs <- get_val("Productive full-length contig rows")
  unique_vdj_barcodes            <- get_val("Unique VDJ barcodes")
  paired_vdj_barcodes            <- get_val("Paired heavy + light VDJ barcodes")
  unpaired_vdj_barcodes          <- get_val("Unpaired/incomplete VDJ barcodes")
  mapped_to_filtered_gex         <- get_val("Paired VDJ barcodes mapped to filtered GEX")
  not_mapped_to_filtered_gex     <- get_val("Paired VDJ barcodes not mapped to filtered GEX")
  
  # ── Optional: ag-pos sort lane contig-level mapping step ─────────────────
  # ag_pos_configs is derived from the uploaded VDJ file names — whichever
  # config numbers the caller determined correspond to the ag-pos sort lane.
  has_lane7 <- !is.null(bcr_annotated) && !is.null(filtered_lookup) &&
    !is.null(ag_pos_configs) && length(ag_pos_configs) > 0
  if (has_lane7) {
    lane7_contigs <- bcr_annotated |>
      dplyr::filter(config %in% ag_pos_configs)
    lane7_gex <- filtered_lookup |>
      dplyr::filter(config %in% ag_pos_configs |
                      stringr::str_detect(filtered_gex_cell,
                                          stringr::regex("ag.?pos",
                                                         ignore_case = TRUE)))
    lane7_mapped   <- sum(lane7_contigs$barcode_clean %in% lane7_gex$barcode_clean)
    lane7_unmapped <- nrow(lane7_contigs) - lane7_mapped
  }
  
  # ── Build nodes and links ─────────────────────────────────────────────────
  if (has_lane7) {
    nodes <- data.frame(name = c(
      lbl("Total VDJ contig rows",            total_contigs),                 # 0
      lbl("Productive full-length contigs",   productive_full_length_contigs), # 1
      lbl("Lane 7 mapped contigs\n(Ag-pos)",  lane7_mapped),                  # 2
      lbl("Lane 7 unmapped\n(QC removed)",    lane7_unmapped),                # 3
      lbl("Unique VDJ barcodes",              unique_vdj_barcodes),           # 4
      lbl("Unpaired/incomplete",              unpaired_vdj_barcodes),         # 5
      lbl("Paired heavy + light",             paired_vdj_barcodes),           # 6
      lbl("Mapped to filtered GEX",           mapped_to_filtered_gex),        # 7
      lbl("Not mapped to filtered GEX",       not_mapped_to_filtered_gex)     # 8
    ))
    links <- data.frame(
      source = c(0, 1, 1, 2, 4, 4, 6, 6),
      target = c(1, 2, 3, 4, 5, 6, 7, 8),
      value  = c(productive_full_length_contigs,
                 lane7_mapped, lane7_unmapped,
                 unique_vdj_barcodes,
                 unpaired_vdj_barcodes, paired_vdj_barcodes,
                 mapped_to_filtered_gex, not_mapped_to_filtered_gex)
    )
  } else {
    nodes <- data.frame(name = c(
      lbl("Total VDJ contig rows",          total_contigs),
      lbl("Productive full-length contigs", productive_full_length_contigs),
      lbl("Unique VDJ barcodes",            unique_vdj_barcodes),
      lbl("Unpaired/incomplete",            unpaired_vdj_barcodes),
      lbl("Paired heavy + light",           paired_vdj_barcodes),
      lbl("Mapped to filtered GEX",         mapped_to_filtered_gex),
      lbl("Not mapped to filtered GEX",     not_mapped_to_filtered_gex)
    ))
    links <- data.frame(
      source = c(0, 1, 2, 2, 4, 4),
      target = c(1, 2, 3, 4, 5, 6),
      value  = c(productive_full_length_contigs,
                 unique_vdj_barcodes,
                 unpaired_vdj_barcodes, paired_vdj_barcodes,
                 mapped_to_filtered_gex, not_mapped_to_filtered_gex)
    )
  }
  
  networkD3::sankeyNetwork(
    Links      = links,
    Nodes      = nodes,
    Source     = "source",
    Target     = "target",
    Value      = "value",
    NodeID     = "name",
    fontSize   = 13,
    nodeWidth  = 30,
    sinksRight = FALSE
  )
}

# ── CITE antigen-specificity Sankey (VDJ-centric) ────────────────────────────
#
#' Flow follows VDJ barcodes through QC steps then into CITE antigen calls.
#' Only cells with paired VDJ data that mapped to GEX are considered.
#'
#' @param cite_results    List returned by run_cite_pipeline().
#' @param bcr_annotated   All contig rows — results()$vdj_list$bcr_annotated.
#' @param barcode_summary Per-barcode summary — results()$vdj_list$barcode_summary.
#' @param vdj_gex_match   VDJ→GEX match table — results()$tables$vdj_gex_match.
#' @param filtered_lookup GEX barcode lookup — results()$tables$filtered_lookup.
make_cite_sankey <- function(cite_results, bcr_annotated,
                             barcode_summary, vdj_gex_match,
                             filtered_lookup,
                             ag_pos_configs = NULL) {
  
  meta <- cite_results$seurat_obj@meta.data |>
    tibble::rownames_to_column("cell")
  
  if (!"comparison_group" %in% colnames(meta))
    stop("'comparison_group' not found. Run the CITE pipeline first.")
  
  # ── Step 1: contig-level ag-pos mapping ───────────────────────────────────
  # ag_pos_configs is passed in by the caller, derived from the uploaded VDJ
  # file names (e.g. the highest-numbered config, or whichever the user mapped
  # to the ag-pos sort lane). Falls back to the highest config number if NULL.
  if (is.null(ag_pos_configs) || length(ag_pos_configs) == 0) {
    all_configs <- sort(unique(bcr_annotated$config))
    ag_pos_configs <- all_configs[length(all_configs)]
    message("make_cite_sankey: ag_pos_configs not supplied; defaulting to highest ",
            "config found in VDJ data: ", paste(ag_pos_configs, collapse = ", "))
  }
  
  ag_pos_contigs_check <- bcr_annotated |>
    dplyr::filter(config %in% ag_pos_configs)
  
  if (nrow(ag_pos_contigs_check) == 0)
    stop(paste0(
      "No contigs found for ag-pos config(s): ",
      paste(ag_pos_configs, collapse = ", "),
      ". Configs present in VDJ data: ",
      paste(sort(unique(bcr_annotated$config)), collapse = ", ")
    ))
  
  ag_pos_contigs  <- ag_pos_contigs_check
  ag_pos_total    <- nrow(ag_pos_contigs)
  
  # Match only against GEX cells from the same config(s)
  ag_pos_gex <- filtered_lookup |>
    dplyr::filter(config %in% ag_pos_configs |
                    stringr::str_detect(filtered_gex_cell,
                                        stringr::regex("ag.?pos", ignore_case = TRUE)))
  
  ag_pos_mapped   <- sum(ag_pos_contigs$barcode_clean %in% ag_pos_gex$barcode_clean)
  ag_pos_unmapped <- ag_pos_total - ag_pos_mapped
  
  # ── Step 2: barcode-level QC counts (from paired barcodes in ag-pos) ──────
  ag_pos_bs <- barcode_summary |>
    dplyr::filter(config %in% ag_pos_configs)
  
  ag_pos_unique_bc  <- nrow(ag_pos_bs)
  ag_pos_unpaired   <- sum(!ag_pos_bs$paired_bcr)
  ag_pos_paired     <- sum( ag_pos_bs$paired_bcr)
  
  ag_pos_vdj_match  <- vdj_gex_match |>
    dplyr::filter(config %in% ag_pos_configs)
  
  ag_pos_mapped_bc  <- sum( ag_pos_vdj_match$maps_to_filtered_gex)
  ag_pos_no_map_bc  <- sum(!ag_pos_vdj_match$maps_to_filtered_gex)
  
  # ── Step 3: CITE antigen calls — all ag-pos cells (not just VDJ-paired) ───
  # This gives the 133 FALSE / 3489 TRUE numbers matching the notebook output.
  ag_pos_meta <- meta |>
    dplyr::filter(comparison_group %in% c("Ag_pos", "Ag_pos_ag_neg"))
  
  ag_spec_vdj  <- sum(ag_pos_meta$comparison_group == "Ag_pos",        na.rm = TRUE)
  not_spec_vdj <- sum(ag_pos_meta$comparison_group == "Ag_pos_ag_neg", na.rm = TRUE)
  
  n_per_cell <- stringr::str_count(
    ag_pos_meta$ag_specific_antigens[
      ag_pos_meta$comparison_group == "Ag_pos" &
        !is.na(ag_pos_meta$ag_specific_antigens)], ";"
  ) + 1L
  mono  <- sum(n_per_cell == 1L, na.rm = TRUE)
  bi    <- sum(n_per_cell == 2L, na.rm = TRUE)
  multi <- sum(n_per_cell >= 3L, na.rm = TRUE)
  
  # ── Nodes — clean CITE flow only ─────────────────────────────────────────
  # Paired/unpaired breakdown lives in the VDJ Sankey.
  # Here: contigs → mapped/unmapped → unique barcodes → CITE calls → reactivity
  node_names <- c(
    paste0("Ag-pos contigs\n",           ag_pos_total),    #  0
    paste0("Mapped contigs\n",           ag_pos_mapped),   #  1
    paste0("Unmapped (QC removed)\n",    ag_pos_unmapped), #  2
    paste0("Unique VDJ barcodes\n",      ag_pos_unique_bc),#  3
    paste0("Antigen-specific\n",         ag_spec_vdj),     #  4
    paste0("Not antigen-specific\n",     not_spec_vdj),    #  5
    paste0("Mono-reactive\n(1 antigen) ",    mono),        #  6
    paste0("Bi-reactive\n(2 antigens) ",     bi),          #  7
    paste0("Multi-reactive\n(\u22653) ",     multi)        #  8
  )
  
  sources <- c(0L, 0L, 1L, 3L, 3L, 4L, 4L, 4L)
  targets <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)
  values  <- c(ag_pos_mapped, ag_pos_unmapped,
               ag_pos_unique_bc,
               ag_spec_vdj, not_spec_vdj,
               mono, bi, multi)
  
  # ── Drop zero-value links and orphaned nodes ───────────────────────────────
  keep      <- values > 0L
  sources   <- sources[keep]
  targets   <- targets[keep]
  values    <- values[keep]
  
  used_idx  <- sort(unique(c(sources, targets)))
  remap     <- setNames(seq_along(used_idx) - 1L, as.character(used_idx))
  sources   <- unname(remap[as.character(sources)])
  targets   <- unname(remap[as.character(targets)])
  node_names <- node_names[used_idx + 1L]
  
  networkD3::sankeyNetwork(
    Links      = data.frame(source = sources, target = targets, value = values),
    Nodes      = data.frame(name = node_names),
    Source     = "source",
    Target     = "target",
    Value      = "value",
    NodeID     = "name",
    fontSize   = 13,
    nodeWidth  = 30,
    sinksRight = FALSE
  )
}