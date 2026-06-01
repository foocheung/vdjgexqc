make_sankey <- function(flow_summary) {
  get_val <- function(metric_name) {
    x <- flow_summary |>
      dplyr::filter(metric == metric_name) |>
      dplyr::pull(value) |>
      as.numeric()
    if (length(x) == 0 || is.na(x)) 0 else x
  }

  total_contigs <- get_val("Total VDJ contig rows")
  productive_full_length_contigs <- get_val("Productive full-length contig rows")
  unique_vdj_barcodes <- get_val("Unique VDJ barcodes")
  paired_vdj_barcodes <- get_val("Paired heavy + light VDJ barcodes")
  unpaired_vdj_barcodes <- get_val("Unpaired/incomplete VDJ barcodes")
  mapped_to_filtered_gex <- get_val("Paired VDJ barcodes mapped to filtered GEX")
  not_mapped_to_filtered_gex <- get_val("Paired VDJ barcodes not mapped to filtered GEX")

  nodes <- data.frame(
    name = c(
      paste0("Total VDJ contig rows\n", total_contigs),
      paste0("Productive full-length contigs\n", productive_full_length_contigs),
      paste0("Unique VDJ barcodes\n", unique_vdj_barcodes),
      paste0("Unpaired/incomplete VDJ barcodes\n", unpaired_vdj_barcodes),
      paste0("Paired heavy + light VDJ barcodes\n", paired_vdj_barcodes),
      paste0("Mapped to filtered GEX\n", mapped_to_filtered_gex),
      paste0("Not mapped to filtered GEX\n", not_mapped_to_filtered_gex)
    )
  )

  links <- data.frame(
    source = c(0, 1, 2, 2, 4, 4),
    target = c(1, 2, 3, 4, 5, 6),
    value = c(
      productive_full_length_contigs,
      unique_vdj_barcodes,
      unpaired_vdj_barcodes,
      paired_vdj_barcodes,
      mapped_to_filtered_gex,
      not_mapped_to_filtered_gex
    )
  )

  networkD3::sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    fontSize = 13,
    nodeWidth = 30,
    sinksRight = FALSE
  )
}
