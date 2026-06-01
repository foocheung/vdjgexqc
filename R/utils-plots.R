make_plots <- function(tables) {
  clone_size_plot <- make_clone_size_plot(tables$clone_sizes)
  v_gene_plot <- make_gene_usage_plot(tables$v_gene_usage, gene_col = "v_gene", title = "Top V genes")
  j_gene_plot <- make_gene_usage_plot(tables$j_gene_usage, gene_col = "j_gene", title = "Top J genes")

  list(
    clone_size_plot = clone_size_plot,
    v_gene_plot = v_gene_plot,
    j_gene_plot = j_gene_plot
  )
}

make_clone_size_plot <- function(clone_sizes, top_n = 25) {
  if (nrow(clone_sizes) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No clonotype data available") +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }

  plot_df <- clone_sizes |>
    dplyr::slice_max(clone_size, n = top_n) |>
    dplyr::mutate(raw_clonotype_id = factor(raw_clonotype_id, levels = rev(raw_clonotype_id)))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = raw_clonotype_id, y = clone_size)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Clonotype",
      y = "Clone size",
      title = paste0("Top ", top_n, " clonotypes by clone size")
    ) +
    ggplot2::theme_minimal(base_size = 13)

  plotly::ggplotly(p)
}

make_gene_usage_plot <- function(gene_usage, gene_col, title, top_n = 25) {
  if (nrow(gene_usage) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No gene usage data available") +
      ggplot2::theme_void()
    return(plotly::ggplotly(p))
  }

  plot_df <- gene_usage |>
    dplyr::group_by(.data[[gene_col]]) |>
    dplyr::summarise(n = sum(n), .groups = "drop") |>
    dplyr::slice_max(n, n = top_n) |>
    dplyr::mutate(gene = factor(.data[[gene_col]], levels = rev(.data[[gene_col]])))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = gene, y = n)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Contigs",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 13)

  plotly::ggplotly(p)
}
