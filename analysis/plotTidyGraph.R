#' @title Visualise a network from a tidygraph object
#' @description Visualise a gene-expression network
#' @details
#' Draws a network graph of genes connected to gene-sets. By default,
#' up-regulated genes will be coloured red, whlst down-regulated genes will be
#' coloured blue. This requires a column labelled 'logFC' to be included in the
#' nodes component of the tidygraph.
#'
#' It is also assumed that an \code{NA} value for logFC will be provided for the
#' gene-set nodes. This is used to generate colours for each node, with edges
#' from each node being drawn in the same colours.
#'
#' Available palettes for nodes can be found using
#' \code{rownames(RColorBrewer::brewer.pal.info)} for the RColorBrewer palettes,
#' or \code{hcl.pals()} for the grDevices palettes
#'
#' @param tg The network as a tidygraph
#' @param layout The layout algorithm to apply
#' @param upCol The colour to label up-regulated genes
#' @param downCol The colour to label down-regulated genes
#' @param palette The name of the palette to use. Can be drawn from those
#' provided in RColorBrewer::brewer.pal, or grDevices::hcl.colors
#' @param paletteSource Choose either RColorBrewer or grDevices palettes
#' @param nodeLabelSize This is an arbitrary value
#' @param gene_repel logical(1). Should the gene labels repel away form the points
#' @param label.padding Set the padding within each label around the text
#'
plotTidyGraph <- function(
  tg, layout = "fr", upCol = "#FF000033", downCol = "#0000FF33",
  palette = "Dark2", paletteSource = c("RColorBrewer", "grDevices"),
  nodeLabelSize = 3, gene_repel = TRUE, label.padding = unit(0.15, "lines")
){

  stopifnot(is(tg, "tbl_graph"))
  stopifnot("logFC" %in% names(as_tibble(activate(tg, nodes))))

  paletteSource <- match.arg(paletteSource)

  # Define the palette
  nNodes <- activate(tg, nodes) %>%
    dplyr::filter(is.na(logFC)) %>%
    as_tibble() %>%
    nrow()



  if (paletteSource == "RColorBrewer"){
    info <- RColorBrewer::brewer.pal.info
    palette <- match.arg(palette, rownames(info))
    maxCols <- info[palette, "maxcolors"]
    n <- min(maxCols, nNodes)
    node_pal <- RColorBrewer::brewer.pal(n, palette)
    node_pal <- colorRampPalette(node_pal)(nNodes)
  }
  if (paletteSource == "grDevices"){
    palette <- match.arg(palette, hcl.pals())
    node_pal <- hcl.colors(nNodes, palette)
  }

  ggraph(tg, layout = layout) +
    # Add the edges, using the source node for colours
    geom_edge_arc(
      aes(color = as.character(from)),
      alpha = 0.5,
      show.legend = FALSE,
      strength = 0.5
    ) +
    # Add the gene-set nodes
    geom_node_point(
      aes(fill = as.character(id)),
      data = . %>%
        dplyr::filter(is.na(logFC)),
      size = 10,
      shape = 21,
      stroke = 0.5,
      show.legend = FALSE
    ) +
    # Labels for gene sets
    geom_node_label(
      aes(label = label, fill = as.character(id)),
      data = . %>%
        dplyr::filter(is.na(logFC)),
      repel = TRUE,
      size = nodeLabelSize,
      force = 0.2,
      label.padding = label.padding,
      alpha = 0.8
    ) +
    # Upregulated genes
    geom_node_point(
      aes(size = abs(logFC)),
      data = . %>% dplyr::filter(logFC > 0),
      fill = upCol,
      shape = 21,
      stroke = 0.5,
      show.legend = FALSE,
      alpha = 0.7
    ) +
    # Down-regulated genes
    geom_node_point(
      aes(size = abs(logFC)),
      data = . %>% dplyr::filter(logFC < 0),
      fill = downCol,
      shape = 21,
      stroke = 0.5,
      show.legend = FALSE,
      alpha = 0.7
    ) +
    # Gene labels
    geom_node_text(
      aes(label = label, size = abs(logFC)),
      data = . %>% dplyr::filter(!is.na(logFC)),
      repel = gene_repel,
      colour = "black"
    ) +
    scale_fill_manual(
      values = node_pal,
      na.value = "grey80"
    ) +
    scale_edge_colour_manual(
      values = node_pal,
      na.value = "grey80"
    ) +
    scale_size_continuous(trans = "sqrt", range = c(1, 3.5)) +
    theme_graph() +
    theme(legend.position = "none")

}
