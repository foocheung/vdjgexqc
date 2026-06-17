#' Application UI
#'
#' @import shiny
app_ui <- function(request) {
  tagList(
    golem::activate_js(),
    golem::favicon(),
    fluidPage(
      titlePanel("VDJ \u2194 GEX QC Dashboard"),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          
          # ── File inputs ──────────────────────────────────────────────────
          h4("Data inputs"),
          fileInput("seurat_filt_rds", "Filtered Seurat object (.rds)", accept = ".rds"),
          fileInput("seurat_merged_rds", "Optional: merged/pre-QC Seurat (.rds)", accept = ".rds"),
          fileInput("vdj_files", "VDJ filtered_contig_annotations.csv",
                    accept = ".csv", multiple = TRUE),
          
          # ── Metadata column selectors ────────────────────────────────────
          tags$hr(),
          h4("Metadata columns"),
          helpText("Populated automatically when a Seurat object is uploaded."),
          selectInput("subject_col",       "Subject column",              choices = NULL),
          selectInput("timepoint_col",     "Timepoint column (optional)", choices = NULL),
          selectInput("condition_col",     "Condition / lane column",     choices = NULL),
          selectInput("celltype_fine_col", "Fine cell-type column",       choices = NULL),
          selectInput("celltype_main_col", "Main cell-type column",       choices = NULL),
          
          # ── Barcode parsing ──────────────────────────────────────────────
          tags$hr(),
          h4("Barcode parsing"),
          textInput("config_regex",  "Config regex",  value = "(?<=multi_config_)[0-9]+"),
          textInput("prefix_regex",  "Prefix to remove", value = "^multi_config_[0-9]+_"),
          textInput("suffix_regex",  "Suffix to remove", value = "-1$"),
          helpText("Default: cells look like multi_config_1_AAAC...-1"),
          
          # ── CITE-seq options ─────────────────────────────────────────────
          tags$hr(),
          h4("CITE-seq / ADT options"),
          checkboxInput("run_cite", "Run CITE-seq antigen scoring", value = FALSE),
          conditionalPanel(
            condition = "input.run_cite == true",
            textInput("adt_assay", "ADT assay name", value = "ADT"),
            
            # Control antigen picker (HSA) — dynamic
            uiOutput("cite_hsa_picker"),
            tags$br(),
            
            # Viral antigen picker — dynamic
            uiOutput("cite_antigen_picker"),
            tags$br(),
            
            tags$hr(),
            h5("Sort lane mapping"),
            helpText("Which values in the condition column represent the antigen-positive and antigen-negative sorts?"),
            selectInput("ag_pos_value", "Antigen-POSITIVE sort value", choices = NULL),
            selectInput("ag_neg_value", "Antigen-NEGATIVE sort value", choices = NULL),
            
            tags$hr(),
            selectInput(
              "cite_threshold", "Specificity threshold",
              choices = c(
                "Loose (0.5 diff, 1.5x)"    = "loose_0.5diff_1.5x",
                "Medium (1.0 diff, 2x)"      = "medium_1diff_2x",
                "Strict (1.5 diff, 3x)"      = "strict_1.5diff_3x",
                "Very strict (2.0 diff, 4x)" = "very_strict_2diff_4x"
              ),
              selected = "medium_1diff_2x"
            )
          ),
          
          tags$hr(),
          actionButton("run", "Run analysis", class = "btn-primary btn-block"),
          
          # ── Downloads ────────────────────────────────────────────────────
          tags$hr(),
          h5("Downloads"),
          downloadButton("download_flow_summary",  "Flow summary"),        br(), br(),
          downloadButton("download_loss_summary",  "Barcode loss"),        br(), br(),
          downloadButton("download_vdj_qc",        "VDJ QC"),              br(), br(),
          downloadButton("download_by_config",     "Per-config"),          br(), br(),
          downloadButton("download_clone_sizes",   "Clone sizes"),         br(), br(),
          downloadButton("download_v_gene",        "V gene usage"),        br(), br(),
          downloadButton("download_j_gene",        "J gene usage"),        br(), br(),
          downloadButton("download_sankey_html",   "Sankey HTML"),         br(), br(),
          conditionalPanel(
            condition = "input.run_cite == true",
            downloadButton("download_cite_scores",     "ADT scores"),      br(), br(),
            downloadButton("download_cite_calls",      "Antigen calls"),   br(), br(),
            downloadButton("download_cite_threshold",  "Threshold sensitivity"), br(), br(),
            downloadButton("download_cite_subject",    "Subject summary"), br(), br(),
            downloadButton("download_cite_crossreact", "Cross-reactivity")
          )
        ),
        
        mainPanel(
          width = 9,
          tabsetPanel(
            
            # ── VDJ / GEX tabs ─────────────────────────────────────────────
            tabPanel("Sankey", br(),
                     networkD3::sankeyNetworkOutput("sankey", height = "650px")),
            tabPanel("Flow summary",   br(), DT::DTOutput("flow_summary")),
            tabPanel("Barcode loss",   br(), DT::DTOutput("barcode_loss_summary")),
            tabPanel("VDJ QC",         br(), DT::DTOutput("vdj_qc_summary")),
            tabPanel("Per config",     br(), DT::DTOutput("by_config")),
            tabPanel("Clone sizes",    br(),
                     plotly::plotlyOutput("clone_size_plot", height = "500px"),
                     br(), DT::DTOutput("clone_sizes")),
            tabPanel("V/J usage",      br(),
                     h4("Top V genes"),
                     plotly::plotlyOutput("v_gene_plot", height = "450px"),
                     br(),
                     h4("Top J genes"),
                     plotly::plotlyOutput("j_gene_plot", height = "450px")),
            tabPanel("Matched cells",  br(), DT::DTOutput("matched_cells")),
            
            # ── CITE-seq / figures tab ──────────────────────────────────────
            tabPanel(
              "CITE-seq / Figures",
              br(),
              conditionalPanel(
                condition = "input.run_cite == false",
                tags$p(tags$em(
                  'Enable "Run CITE-seq antigen scoring" in the sidebar and click Run.'
                ))
              ),
              conditionalPanel(
                condition = "input.run_cite == true",
                tabsetPanel(
                  
                  # Fig 1 ─────────────────────────────────────────────────
                  tabPanel(
                    "Fig 1 \u2014 Subject cell counts",
                    br(),
                    tags$p("Cells per subject after QC (replicates slide 1)."),
                    DT::DTOutput("fig_subject_table")
                  ),
                  
                  # Fig 3 ─────────────────────────────────────────────────
                  tabPanel(
                    "Fig 3 \u2014 ADT dot plot",
                    br(),
                    tags$p("ADT average expression \u00d7 % expressed per cell type (replicates Q1d dot plot)."),
                    plotly::plotlyOutput("fig_adt_dotplot", height = "600px")
                  ),
                  
                  # Fig 4 ─────────────────────────────────────────────────
                  tabPanel(
                    "Fig 4 \u2014 B-cell abundance",
                    br(),
                    tags$p("B-cell subpopulation abundance as % of all PBMCs per subject."),
                    plotly::plotlyOutput("fig_bcell_abundance", height = "500px")
                  ),
                  
                  # Fig 5 ─────────────────────────────────────────────────
                  tabPanel(
                    "Fig 5 \u2014 Cell-type heatmap",
                    br(),
                    tags$p("Z-scored cell-type frequencies per subject (clustered)."),
                    plotly::plotlyOutput("fig_celltype_heatmap", height = "600px")
                  ),
                  
                  # Fig 7 ─────────────────────────────────────────────────
                  tabPanel(
                    "Fig 7 \u2014 Antigen vs HSA",
                    br(),
                    tags$p("Per-cell antigen signal vs control antigen (HSA) faceted by subject and antigen. Red = above specificity threshold."),
                    plotly::plotlyOutput("fig_ag_vs_hsa", height = "700px")
                  ),
                  
                  # Fig 8 ─────────────────────────────────────────────────
                  tabPanel(
                    "Fig 8 \u2014 B-cells by group",
                    br(),
                    tags$p("B-cell subtype % within Ag-neg sort vs antigen-positive cells, per subject."),
                    plotly::plotlyOutput("fig_bcell_by_group", height = "600px")
                  ),
                  
                  # Fig 9 ─────────────────────────────────────────────────
                  tabPanel(
                    "Fig 9 \u2014 IGH clonotypes",
                    br(),
                    tags$p("Top IGH clonotypes in Ag-neg vs antigen-positive lanes for a chosen subject."),
                    uiOutput("clonotype_subject_ui"),
                    br(),
                    plotly::plotlyOutput("fig_clonotype_lollipop", height = "700px")
                  ),
                  
                  # CITE QC sub-tabs ───────────────────────────────────────
                  tabPanel(
                    "CITE QC",
                    br(),
                    h4("Detected ADT features"),
                    verbatimTextOutput("cite_adt_features"),
                    tags$hr(),
                    h4("Threshold sensitivity"),
                    DT::DTOutput("cite_threshold_sensitivity"),
                    tags$hr(),
                    h4("Subject summary"),
                    DT::DTOutput("cite_subject_summary"),
                    tags$hr(),
                    h4("Cross-reactivity"),
                    DT::DTOutput("cite_cross_react"),
                    tags$hr(),
                    h4("ADT score distributions"),
                    plotly::plotlyOutput("cite_score_plot", height = "450px"),
                    tags$hr(),
                    h4("Antigen calls (chosen threshold)"),
                    DT::DTOutput("cite_calls_table")
                  )
                )
              )
            ),
            
            tabPanel("Notes", br(),
                     tags$p("Each VDJ row is a contig/receptor chain."),
                     tags$p("A paired BCR barcode has at least one IGH and one IGK/IGL."),
                     tags$p("Mapped to GEX = cleaned VDJ barcode + config found in the filtered Seurat object."),
                     tags$p("CITE-seq: antigen-specific calls require signal > control antigen by both the absolute difference and fold-ratio thresholds selected in the sidebar."),
                     tags$p("Fig 9 requires VDJ data with a cdr3, cdr3_aa, or junction_aa column."))
          )
        )
      )
    )
  )
}