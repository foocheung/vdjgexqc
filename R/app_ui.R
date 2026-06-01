#' Application UI
#'
#' @import shiny
app_ui <- function(request) {
  tagList(
    golem::activate_js(),
    golem::favicon(),
    fluidPage(
      titlePanel("VDJ ↔ GEX QC Dashboard"),
      sidebarLayout(
        sidebarPanel(
          fileInput(
            "seurat_filt_rds",
            "Upload filtered/final Seurat object (.rds)",
            accept = c(".rds")
          ),
          fileInput(
            "seurat_merged_rds",
            "Optional: Upload merged/pre-QC Seurat object (.rds)",
            accept = c(".rds")
          ),
          fileInput(
            "vdj_files",
            "Upload VDJ filtered_contig_annotations.csv files",
            accept = c(".csv"),
            multiple = TRUE
          ),

          tags$hr(),
          h4("Barcode parsing options"),

          textInput(
            "config_regex",
            "Regex to extract config from Seurat cell names",
            value = "(?<=multi_config_)[0-9]+"
          ),
          textInput(
            "prefix_regex",
            "Regex prefix to remove from Seurat cell names",
            value = "^multi_config_[0-9]+_"
          ),
          textInput(
            "suffix_regex",
            "Regex suffix to remove from Seurat cell names",
            value = "-1$"
          ),

          helpText(
            "Default assumes Seurat cells look like multi_config_1_AAAC...-1 and VDJ files look like multi_config_1_filtered_contig_annotations.csv."
          ),

          actionButton("run", "Run VDJ ↔ GEX QC", class = "btn-primary"),

          tags$hr(),
          downloadButton("download_flow_summary", "Download flow summary CSV"),
          br(), br(),
          downloadButton("download_loss_summary", "Download barcode loss CSV"),
          br(), br(),
          downloadButton("download_vdj_qc", "Download VDJ QC CSV"),
          br(), br(),
          downloadButton("download_by_config", "Download per-config CSV"),
          br(), br(),
          downloadButton("download_clone_sizes", "Download clone sizes CSV"),
          br(), br(),
          downloadButton("download_v_gene", "Download V gene usage CSV"),
          br(), br(),
          downloadButton("download_j_gene", "Download J gene usage CSV"),
          br(), br(),
          downloadButton("download_sankey_html", "Download Sankey HTML")
        ),

        mainPanel(
          tabsetPanel(
            tabPanel(
              "Sankey",
              br(),
              networkD3::sankeyNetworkOutput("sankey", height = "650px")
            ),
            tabPanel(
              "Flow summary",
              br(),
              DT::DTOutput("flow_summary")
            ),
            tabPanel(
              "Barcode loss",
              br(),
              DT::DTOutput("barcode_loss_summary")
            ),
            tabPanel(
              "VDJ QC summary",
              br(),
              DT::DTOutput("vdj_qc_summary")
            ),
            tabPanel(
              "Per config",
              br(),
              DT::DTOutput("by_config")
            ),
            tabPanel(
              "Clone sizes",
              br(),
              plotly::plotlyOutput("clone_size_plot", height = "500px"),
              br(),
              DT::DTOutput("clone_sizes")
            ),
            tabPanel(
              "V/J usage",
              br(),
              h4("Top V genes"),
              plotly::plotlyOutput("v_gene_plot", height = "450px"),
              br(),
              h4("Top J genes"),
              plotly::plotlyOutput("j_gene_plot", height = "450px")
            ),
            tabPanel(
              "Matched cells",
              br(),
              DT::DTOutput("matched_cells")
            ),
            tabPanel(
              "Notes",
              br(),
              tags$p("Each VDJ row is a contig/receptor chain."),
              tags$p("Barcode means a cell barcode. A paired BCR barcode is a cell barcode with at least one IGH and at least one IGK or IGL."),
              tags$p("Mapped to GEX means the cleaned VDJ barcode and config were found in the uploaded filtered Seurat object."),
              tags$p("If a merged/pre-QC Seurat object is uploaded, the Barcode loss tab tracks paired VDJ barcodes through merged and filtered objects.")
            )
          )
        )
      )
    )
  )
}
