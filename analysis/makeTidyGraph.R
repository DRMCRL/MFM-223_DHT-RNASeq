#' @title Form a tidygraph object
#' @description Form a tidygraph from a list of genesets
#' @details
#' Taking a list of genesets as the primary input form a network. Names of each
#' list element are assumed to be the names of genesets, with the values within
#' each element assumed to be gene identifiers. Choosing actual gene names may
#' be the most useful option.
#'
#' If provided, the topTable will be joined onto the nodes allowing all columns
#' to be used for modifying the final plot. Again, the most useful columns may
#' be either logFC or the PValue.
#'
#' By default, all nodes (i.e. gene-sets and gene names) will be contained in a
#' column named 'label'. Setting the correct column from the topTable to join
#' on can be performed using the .by argument
#' @param geneSets A list of gene-sets as described in the details section
#' @param topTable An optional topTable as output by a function such as
#' limma::topTable
#' @param .by If providing the
makeTidyGraph <- function(geneSets, topTable, .by = c("label" = "gene_name")){

  ## Create a node list
  nodes <- tibble(label = c(names(geneSets), unlist(geneSets))) %>%
    distinct(label) %>%
    rowid_to_column("id")

  # Then create an edge list
  edges <- plyr::ldply(geneSets, tibble) %>%
    set_colnames(c("pathway", "gene")) %>%
    left_join(nodes, by = c("pathway" = "label")) %>%
    dplyr::rename(from = id) %>%
    left_join(nodes, by = c("gene" = "label")) %>%
    dplyr::rename(to = id) %>%
    dplyr::select(from, to)

  g <- tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = FALSE
  )
  if (!missing(topTable)){
    g <- activate(g, nodes) %>%
      left_join(topTable, by = .by)
  }
  g

}
