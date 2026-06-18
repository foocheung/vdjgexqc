#' Application server
#'
#' @import shiny
app_server <- function(input, output, session) {
  
  # ── Load Seurat on upload ───────────────────────────────────────────────────
  seurat_filt_obj <- reactive({
    req(input$seurat_filt_rds)
    readRDS(input$seurat_filt_rds$datapath)
  })
  
  # ── Populate metadata column pickers once Seurat is loaded ──────────────────
  observe({
    req(seurat_filt_obj())
    meta_cols <- colnames(seurat_filt_obj()@meta.data)
    
    updateSelectInput(session, "subject_col",
                      choices  = meta_cols,
                      selected = intersect(c("subject", "orig.ident", meta_cols[1]), meta_cols)[1])
    
    updateSelectInput(session, "timepoint_col",
                      choices  = c("(none)" = "", meta_cols),
                      selected = intersect(c("timepoint", ""), meta_cols)[1])
    
    fine_guess <- intersect(c("monaco_fine", "celltype_fine", meta_cols[1]), meta_cols)[1]
    main_guess <- intersect(c("monaco_main", "celltype_main", meta_cols[1]), meta_cols)[1]
    updateSelectInput(session, "celltype_fine_col", choices = meta_cols, selected = fine_guess)
    updateSelectInput(session, "celltype_main_col", choices = meta_cols, selected = main_guess)
    
    cond_guess <- intersect(c("condition", "lane", "sort", meta_cols[1]), meta_cols)[1]
    updateSelectInput(session, "condition_col", choices = meta_cols, selected = cond_guess)
  })
  
  # ── When condition column changes, repopulate the value pickers ──────────────
  observe({
    req(seurat_filt_obj(), input$condition_col)
    vals <- sort(unique(as.character(
      seurat_filt_obj()@meta.data[[input$condition_col]]
    )))
    vals <- vals[!is.na(vals)]
    
    pos_guess <- intersect(c("Antigen positive", "ag_pos", "positive", vals[1]), vals)[1]
    neg_guess <- intersect(c("Antigen negative", "ag_neg", "negative", vals[length(vals)]), vals)[1]
    
    updateSelectInput(session, "ag_pos_value", choices = vals, selected = pos_guess)
    updateSelectInput(session, "ag_neg_value", choices = vals, selected = neg_guess)
  })
  
  # ── Detect ADT features on upload ───────────────────────────────────────────
  adt_info <- reactive({
    req(seurat_filt_obj(), input$run_cite, input$adt_assay)
    detect_adt_features(seurat_filt_obj(), assay = input$adt_assay)
  })
  
  # ── Dynamic antigen picker (viral antigens) ──────────────────────────────────
  output$cite_antigen_picker <- renderUI({
    info <- adt_info()
    if (is.null(info) || length(info$antigen_features) == 0)
      return(tags$p(tags$em("No non-HSA ADT features found. Check assay name.")))
    
    tagList(
      tags$label(paste0("Viral antigens to score (",
                        length(info$antigen_features), " detected):")),
      tags$div(
        actionLink("cite_select_all",   "Select all"),
        tags$span(" | "),
        actionLink("cite_deselect_all", "Deselect all")
      ),
      checkboxGroupInput("cite_selected_antigens", label = NULL,
                         choices  = info$antigen_features,
                         selected = {
                           dv <- grep("DV[1-4]|DV_[1-4]|dengue.?[1-4]|serotype.?[1-4]",
                                      info$antigen_features,
                                      value = TRUE, ignore.case = TRUE)
                           if (length(dv) > 0) dv else info$antigen_features
                         })
    )
  })
  
  # ── Dynamic HSA / control antigen picker ────────────────────────────────────
  output$cite_hsa_picker <- renderUI({
    info <- adt_info()
    if (is.null(info)) return(NULL)
    
    tagList(
      tags$label("Control antigen(s) — HSA or equivalent:"),
      checkboxGroupInput("cite_hsa_features", label = NULL,
                         choices  = info$all_features,
                         selected = info$hsa_features)
    )
  })
  
  observeEvent(input$cite_select_all, {
    info <- adt_info(); req(info)
    updateCheckboxGroupInput(session, "cite_selected_antigens",
                             selected = info$antigen_features)
  })
  observeEvent(input$cite_deselect_all, {
    updateCheckboxGroupInput(session, "cite_selected_antigens", selected = character(0))
  })
  
  # Show detected features summary in CITE tab
  output$cite_adt_features <- renderText({
    info <- adt_info()
    if (is.null(info)) return("Upload a Seurat object to inspect ADT features.")
    paste0(
      "All ADT features (", length(info$all_features), "): ",
      paste(info$all_features, collapse = ", "),
      "\nControl/HSA feature(s): ", paste(info$hsa_features, collapse = ", "),
      "\nAvailable antigens (", length(info$antigen_features), "): ",
      paste(info$antigen_features, collapse = ", ")
    )
  })
  
  # ── Core VDJ / GEX pipeline ─────────────────────────────────────────────────
  results <- eventReactive(input$run, {
    req(input$seurat_filt_rds, input$vdj_files)
    
    withProgress(message = "Reading Seurat and VDJ files...", value = 0, {
      incProgress(0.15)
      seurat_filt   <- seurat_filt_obj()
      seurat_merged <- NULL
      if (!is.null(input$seurat_merged_rds))
        seurat_merged <- readRDS(input$seurat_merged_rds$datapath)
      
      incProgress(0.35)
      vdj_list <- build_vdj_summary(input$vdj_files)
      
      incProgress(0.60)
      tables <- make_flow_tables(
        vdj_list      = vdj_list,
        seurat_filt   = seurat_filt,
        seurat_merged = seurat_merged,
        config_regex  = input$config_regex,
        prefix_regex  = input$prefix_regex,
        suffix_regex  = input$suffix_regex
      )
      
      # ── Derive ag-pos config(s) from uploaded VDJ file names ────────────────
      # Parse config numbers from every uploaded filename, then identify which
      # correspond to the ag-pos sort lane.  The app UI lets the user specify
      # the ag-pos sort value; we match that against the condition column to
      # find the right config(s).  As a reliable fallback we also expose the
      # parsed config numbers so the server can pass them directly.
      all_file_configs <- purrr::map_chr(input$vdj_files$name, parse_config_from_file)
      ag_pos_configs <- derive_ag_pos_configs(
        vdj_list         = vdj_list,
        file_names       = input$vdj_files$name,
        all_file_configs = all_file_configs,
        ag_pos_value     = input$ag_pos_value,
        seurat_filt      = seurat_filt,
        condition_col    = input$condition_col
      )
      
      incProgress(0.85)
      sankey <- make_sankey(
        flow_summary    = tables$flow_summary,
        bcr_annotated   = vdj_list$bcr_annotated,
        filtered_lookup = tables$filtered_lookup,
        ag_pos_configs  = ag_pos_configs
      )
      plots  <- make_plots(tables)
      incProgress(1)
      
      list(seurat_filt   = seurat_filt,
           seurat_merged  = seurat_merged,
           vdj_list      = vdj_list,
           tables        = tables,
           sankey        = sankey,
           plots         = plots,
           ag_pos_configs = ag_pos_configs)
    })
  })
  
  # ── CITE-seq pipeline ────────────────────────────────────────────────────────
  cite_results <- eventReactive(input$run, {
    req(results())
    if (!isTRUE(input$run_cite)) return(NULL)
    
    selected  <- input$cite_selected_antigens
    hsa_feats <- input$cite_hsa_features
    
    if (is.null(selected) || length(selected) == 0) {
      showNotification("No viral antigens selected.", type = "warning", duration = 8)
      return(NULL)
    }
    if (is.null(hsa_feats) || length(hsa_feats) == 0) {
      showNotification("No control antigen selected.", type = "warning", duration = 8)
      return(NULL)
    }
    
    withProgress(message = "Running CITE-seq antigen scoring...", value = 0, {
      incProgress(0.3)
      out <- tryCatch(
        run_cite_pipeline(
          seurat_obj            = results()$seurat_filt,
          selected_antigens     = selected,
          hsa_features_override = hsa_feats,
          chosen_threshold      = input$cite_threshold,
          adt_assay             = input$adt_assay,
          condition_col         = input$condition_col,
          ag_pos_value          = input$ag_pos_value,
          ag_neg_value          = input$ag_neg_value
        ),
        error = function(e) {
          showNotification(paste("CITE error:", conditionMessage(e)),
                           type = "error", duration = 10)
          NULL
        }
      )
      incProgress(1)
      out
    })
  })
  
  # Helper: subject column from user selection
  subj_col <- reactive({ input$subject_col })
  tp_col   <- reactive({
    v <- input$timepoint_col
    if (is.null(v) || v == "") NULL else v
  })
  
  # ── VDJ / GEX Sankey ────────────────────────────────────────────────────────
  output$sankey <- networkD3::renderSankeyNetwork({ req(results()); results()$sankey })
  
  # ── CITE antigen-filtering Sankey ────────────────────────────────────────────
  output$cite_sankey <- networkD3::renderSankeyNetwork({
    req(cite_results(), results())
    tryCatch(
      make_cite_sankey(
        cite_results    = cite_results(),
        bcr_annotated   = results()$vdj_list$bcr_annotated,
        barcode_summary = results()$vdj_list$barcode_summary,
        vdj_gex_match   = results()$tables$vdj_gex_match,
        filtered_lookup = results()$tables$filtered_lookup,
        ag_pos_configs  = results()$ag_pos_configs
      ),
      error = function(e) {
        showNotification(paste("CITE Sankey error:", conditionMessage(e)),
                         type = "error", duration = 10)
        NULL
      }
    )
  })
  
  output$flow_summary <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$flow_summary, rownames = FALSE,
                  options = list(pageLength = 20))
  })
  output$barcode_loss_summary <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$barcode_loss_summary, rownames = FALSE,
                  options = list(pageLength = 20))
  })
  output$vdj_qc_summary <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$vdj_qc_summary, rownames = FALSE,
                  options = list(pageLength = 20))
  })
  output$by_config <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$by_config, rownames = FALSE,
                  options = list(pageLength = 20))
  })
  output$clone_sizes <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$clone_sizes, rownames = FALSE,
                  options = list(pageLength = 25))
  })
  output$matched_cells <- DT::renderDT({
    req(results())
    DT::datatable(
      results()$tables$vdj_gex_match |>
        dplyr::select(config, barcode, barcode_clean,
                      has_heavy, has_light, n_contigs, n_heavy, n_light,
                      paired_bcr, maps_to_filtered_gex, filtered_gex_cell,
                      maps_to_merged_gex, merged_gex_cell),
      rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE)
    )
  })
  output$clone_size_plot <- plotly::renderPlotly({ req(results()); results()$plots$clone_size_plot })
  output$v_gene_plot     <- plotly::renderPlotly({ req(results()); results()$plots$v_gene_plot })
  output$j_gene_plot     <- plotly::renderPlotly({ req(results()); results()$plots$j_gene_plot })
  
  # ── CITE tab: core outputs ────────────────────────────────────────────────────
  output$cite_threshold_sensitivity <- DT::renderDT({
    req(cite_results())
    DT::datatable(cite_results()$threshold_sensitivity, rownames = FALSE,
                  options = list(pageLength = 20))
  })
  output$cite_subject_summary <- DT::renderDT({
    req(cite_results())
    tbl <- cite_results()$ag_subject_summary
    if (is.null(tbl)) return(DT::datatable(data.frame(message = "No subject/condition metadata.")))
    DT::datatable(tbl, rownames = FALSE, options = list(pageLength = 20))
  })
  output$cite_cross_react <- DT::renderDT({
    req(cite_results())
    tbl <- cite_results()$cross_react
    if (is.null(tbl)) return(DT::datatable(data.frame(message = "No cross-reactivity data.")))
    DT::datatable(tbl, rownames = FALSE, options = list(pageLength = 20))
  })
  output$cite_calls_table <- DT::renderDT({
    req(cite_results())
    DT::datatable(
      cite_results()$ag_calls_wide |>
        dplyr::filter(threshold_label == cite_results()$chosen_threshold),
      rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE)
    )
  })
  output$cite_score_plot <- plotly::renderPlotly({
    req(cite_results())
    scores   <- cite_results()$ag_scores
    antigens <- cite_results()$selected_antigens
    long <- scores |>
      tidyr::pivot_longer(dplyr::all_of(c(antigens, "HSA")),
                          names_to = "antigen", values_to = "score") |>
      dplyr::filter(!is.na(score))
    p <- ggplot2::ggplot(long, ggplot2::aes(x = score, colour = antigen, fill = antigen)) +
      ggplot2::geom_density(alpha = 0.25) +
      ggplot2::facet_wrap(~antigen, scales = "free_y") +
      ggplot2::labs(x = "Normalised ADT signal", y = "Density") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "none")
    plotly::ggplotly(p)
  })
  
  # ── Figure outputs ────────────────────────────────────────────────────────────
  output$fig_subject_table <- DT::renderDT({
    req(results(), subj_col())
    tbl <- tryCatch(
      make_subject_cell_table(results()$seurat_filt, subj_col(), tp_col()),
      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
    )
    req(tbl)
    DT::datatable(tbl, rownames = FALSE, options = list(pageLength = 20))
  })
  
  output$fig_adt_dotplot <- plotly::renderPlotly({
    req(cite_results(), input$celltype_main_col)
    tryCatch(
      make_adt_dotplot(
        seurat_obj       = results()$seurat_filt,
        antigen_features = cite_results()$selected_antigens,
        celltype_col     = input$celltype_main_col,
        adt_assay        = input$adt_assay
      ),
      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
    )
  })
  
  output$fig_bcell_abundance <- plotly::renderPlotly({
    req(results(), subj_col(), input$celltype_fine_col)
    tryCatch(
      make_bcell_abundance_plot(
        seurat_obj    = results()$seurat_filt,
        subject_col   = subj_col(),
        timepoint_col = tp_col(),
        bcell_col     = input$celltype_fine_col
      ),
      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
    )
  })
  
  output$fig_celltype_heatmap <- plotly::renderPlotly({
    req(results(), subj_col(), input$celltype_fine_col)
    tryCatch(
      make_celltype_heatmap(
        seurat_obj   = results()$seurat_filt,
        subject_col  = subj_col(),
        celltype_col = input$celltype_fine_col
      ),
      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
    )
  })
  
  output$fig_ag_vs_hsa <- plotly::renderPlotly({
    req(cite_results(), subj_col())
    tryCatch(
      make_ag_vs_hsa_plot(
        ag_scores        = cite_results()$ag_scores,
        antigen_cols     = cite_results()$selected_antigens,
        subject_col      = subj_col(),
        ag_long          = cite_results()$ag_long,
        chosen_threshold = cite_results()$chosen_threshold
      ),
      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
    )
  })
  
  output$fig_bcell_by_group <- plotly::renderPlotly({
    req(cite_results(), subj_col(), input$celltype_fine_col)
    tryCatch(
      make_bcell_by_group_plot(
        seurat_obj  = cite_results()$seurat_obj,
        subject_col = subj_col(),
        bcell_col   = "bcell_type"
      ),
      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
    )
  })
  
  output$clonotype_subject_ui <- renderUI({
    req(results())
    meta     <- results()$seurat_filt@meta.data
    subj_col <- subj_col()
    subjects <- if (subj_col %in% colnames(meta)) sort(unique(meta[[subj_col]])) else character(0)
    selectInput("clonotype_subject", "Subject for clonotype plot:", choices = subjects)
  })
  
  output$fig_clonotype_lollipop <- plotly::renderPlotly({
    req(cite_results(), results(), input$clonotype_subject)
    tryCatch(
      make_clonotype_lollipop(
        bcr_productive = results()$vdj_list$bcr_productive,
        seurat_meta    = cite_results()$seurat_obj@meta.data,
        subject_id     = input$clonotype_subject,
        subject_col    = subj_col(),
        ag_pos_label   = "Ag_pos",
        ag_neg_label   = "Ag_neg_sort"
      ),
      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }
    )
  })
  
  # ── Downloads ─────────────────────────────────────────────────────────────────
  output$download_flow_summary <- downloadHandler(
    filename = function() "VDJ_to_GEX_flow_summary.csv",
    content  = function(f) { req(results()); readr::write_csv(results()$tables$flow_summary, f) }
  )
  output$download_loss_summary <- downloadHandler(
    filename = function() "VDJ_to_GEX_barcode_loss_summary.csv",
    content  = function(f) { req(results()); readr::write_csv(results()$tables$barcode_loss_summary, f) }
  )
  output$download_vdj_qc <- downloadHandler(
    filename = function() "VDJ_QC_summary.csv",
    content  = function(f) { req(results()); readr::write_csv(results()$tables$vdj_qc_summary, f) }
  )
  output$download_by_config <- downloadHandler(
    filename = function() "VDJ_to_GEX_by_config.csv",
    content  = function(f) { req(results()); readr::write_csv(results()$tables$by_config, f) }
  )
  output$download_clone_sizes <- downloadHandler(
    filename = function() "VDJ_clone_sizes.csv",
    content  = function(f) { req(results()); readr::write_csv(results()$tables$clone_sizes, f) }
  )
  output$download_v_gene <- downloadHandler(
    filename = function() "VDJ_V_gene_usage.csv",
    content  = function(f) { req(results()); readr::write_csv(results()$tables$v_gene_usage, f) }
  )
  output$download_j_gene <- downloadHandler(
    filename = function() "VDJ_J gene_usage.csv",
    content  = function(f) { req(results()); readr::write_csv(results()$tables$j_gene_usage, f) }
  )
  output$download_sankey_html <- downloadHandler(
    filename = function() "VDJ_to_GEX_sankey.html",
    content  = function(f) {
      req(results())
      htmlwidgets::saveWidget(results()$sankey, f, selfcontained = TRUE)
    }
  )
  output$download_cite_scores <- downloadHandler(
    filename = function() "CITE_ADT_scores.csv",
    content  = function(f) { req(cite_results()); readr::write_csv(cite_results()$ag_scores, f) }
  )
  output$download_cite_calls <- downloadHandler(
    filename = function() "CITE_antigen_calls.csv",
    content  = function(f) { req(cite_results()); readr::write_csv(cite_results()$ag_calls_wide, f) }
  )
  output$download_cite_threshold <- downloadHandler(
    filename = function() "CITE_threshold_sensitivity.csv",
    content  = function(f) { req(cite_results()); readr::write_csv(cite_results()$threshold_sensitivity, f) }
  )
  output$download_cite_subject <- downloadHandler(
    filename = function() "CITE_subject_summary.csv",
    content  = function(f) {
      req(cite_results())
      tbl <- cite_results()$ag_subject_summary
      if (!is.null(tbl)) readr::write_csv(tbl, f)
    }
  )
  output$download_annotated_seurat <- downloadHandler(
    filename = function() "seurat_with_DV_calls.rds",
    content  = function(f) {
      req(cite_results())
      saveRDS(cite_results()$seurat_obj, f)
    }
  )
  output$download_cite_crossreact <- downloadHandler(
    filename = function() "CITE_cross_reactivity.csv",
    content  = function(f) {
      req(cite_results())
      tbl <- cite_results()$cross_react
      if (!is.null(tbl)) readr::write_csv(tbl, f)
    }
  )
}