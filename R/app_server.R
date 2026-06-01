#' Application server
#'
#' @import shiny
app_server <- function(input, output, session) {

  results <- eventReactive(input$run, {
    req(input$seurat_filt_rds)
    req(input$vdj_files)

    withProgress(message = "Reading Seurat and VDJ files...", value = 0, {
      incProgress(0.15)

      seurat_filt <- readRDS(input$seurat_filt_rds$datapath)

      seurat_merged <- NULL
      if (!is.null(input$seurat_merged_rds)) {
        seurat_merged <- readRDS(input$seurat_merged_rds$datapath)
      }

      incProgress(0.35)

      vdj_list <- build_vdj_summary(input$vdj_files)

      incProgress(0.60)

      tables <- make_flow_tables(
        vdj_list = vdj_list,
        seurat_filt = seurat_filt,
        seurat_merged = seurat_merged,
        config_regex = input$config_regex,
        prefix_regex = input$prefix_regex,
        suffix_regex = input$suffix_regex
      )

      incProgress(0.85)

      sankey <- make_sankey(tables$flow_summary)
      plots <- make_plots(tables)

      incProgress(1)

      list(
        seurat_filt = seurat_filt,
        seurat_merged = seurat_merged,
        vdj_list = vdj_list,
        tables = tables,
        sankey = sankey,
        plots = plots
      )
    })
  })

  output$sankey <- networkD3::renderSankeyNetwork({
    req(results())
    results()$sankey
  })

  output$flow_summary <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$flow_summary, rownames = FALSE, options = list(pageLength = 20))
  })

  output$barcode_loss_summary <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$barcode_loss_summary, rownames = FALSE, options = list(pageLength = 20))
  })

  output$vdj_qc_summary <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$vdj_qc_summary, rownames = FALSE, options = list(pageLength = 20))
  })

  output$by_config <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$by_config, rownames = FALSE, options = list(pageLength = 20))
  })

  output$clone_sizes <- DT::renderDT({
    req(results())
    DT::datatable(results()$tables$clone_sizes, rownames = FALSE, options = list(pageLength = 25))
  })

  output$matched_cells <- DT::renderDT({
    req(results())
    DT::datatable(
      results()$tables$vdj_gex_match |>
        dplyr::select(
          config,
          barcode,
          barcode_clean,
          has_heavy,
          has_light,
          n_contigs,
          n_heavy,
          n_light,
          paired_bcr,
          maps_to_filtered_gex,
          filtered_gex_cell,
          maps_to_merged_gex,
          merged_gex_cell
        ),
      rownames = FALSE,
      options = list(pageLength = 25, scrollX = TRUE)
    )
  })

  output$clone_size_plot <- plotly::renderPlotly({
    req(results())
    results()$plots$clone_size_plot
  })

  output$v_gene_plot <- plotly::renderPlotly({
    req(results())
    results()$plots$v_gene_plot
  })

  output$j_gene_plot <- plotly::renderPlotly({
    req(results())
    results()$plots$j_gene_plot
  })

  output$download_flow_summary <- downloadHandler(
    filename = function() "VDJ_to_GEX_flow_summary.csv",
    content = function(file) {
      req(results())
      readr::write_csv(results()$tables$flow_summary, file)
    }
  )

  output$download_loss_summary <- downloadHandler(
    filename = function() "VDJ_to_GEX_barcode_loss_summary.csv",
    content = function(file) {
      req(results())
      readr::write_csv(results()$tables$barcode_loss_summary, file)
    }
  )

  output$download_vdj_qc <- downloadHandler(
    filename = function() "VDJ_QC_summary.csv",
    content = function(file) {
      req(results())
      readr::write_csv(results()$tables$vdj_qc_summary, file)
    }
  )

  output$download_by_config <- downloadHandler(
    filename = function() "VDJ_to_GEX_by_config.csv",
    content = function(file) {
      req(results())
      readr::write_csv(results()$tables$by_config, file)
    }
  )

  output$download_clone_sizes <- downloadHandler(
    filename = function() "VDJ_clone_sizes.csv",
    content = function(file) {
      req(results())
      readr::write_csv(results()$tables$clone_sizes, file)
    }
  )

  output$download_v_gene <- downloadHandler(
    filename = function() "VDJ_V_gene_usage.csv",
    content = function(file) {
      req(results())
      readr::write_csv(results()$tables$v_gene_usage, file)
    }
  )

  output$download_j_gene <- downloadHandler(
    filename = function() "VDJ_J_gene_usage.csv",
    content = function(file) {
      req(results())
      readr::write_csv(results()$tables$j_gene_usage, file)
    }
  )

  output$download_sankey_html <- downloadHandler(
    filename = function() "VDJ_to_GEX_sankey.html",
    content = function(file) {
      req(results())
      htmlwidgets::saveWidget(results()$sankey, file = file, selfcontained = TRUE)
    }
  )
}
