#' @include getters_and_setters.R
NULL

#' Run Force-Directed Graph (FDG) on Affinity Matrix
#'
#' This function constructs a force-directed graph to visualize the relationships between genes
#' in selected modules. It uses the igraph package to create an undirected weighted graph and
#' determines the layout using the Fruchterman-Reingold algorithm.
#'
#' @param object A Seurat object.
#' @param module.source A character string indicating the source of the module data within the
#' Seurat object. Default is "GeneCooc".
#' @param module.size An integer indicating the minimum number of genes a module must contain
#' to be included in the analysis. Default is 10.
#' @param weight.cutoff A numeric threshold for the values in the affinity matrix. Values with
#' below this cutoff will be discarded. Default is 0.1.
#' @param on.module An optional character string specifying a particular module to analyze.
#' If NULL, all modules are considered. Default is NULL.
#'
#' @return A modified Seurat object that includes the force-directed graph layout coordinates
#' and adjacency data in its `misc` slot under a list named by `module.source`.
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#'
#' @export
RunModuleFDG <- function(
    object,
    module.source="GeneCooc",
    module.size = 10,
    weight.cutoff = 0.1,
    on.module = NULL
){
  module.list <- FetchModuleList(object, module.source = module.source, module.size = module.size)
  module.DF <- FetchModuleDF(object, module.source = module.source)
  A <- FetchAffinityMatrix(object, module.source = module.source)
  if (is.null(on.module)) {
    genes.use <- unique(base::Reduce(c, module.list))
  } else {
    modules <- names(module.list)[startsWith(names(module.list), prefix = on.module)]
    genes.use <- unique(base::Reduce(c, module.list[modules[-1]]))
  }
  A <- A[genes.use, genes.use]
  A[A < weight.cutoff] <- 0
  g <- igraph::graph_from_adjacency_matrix(adjmatrix = A,
                                           mode = "undirected",
                                           weighted = TRUE,
                                           add.colnames = TRUE)
  fr.layout <- igraph::layout_with_fr(g)
  fr.layout <- as.data.frame(fr.layout)
  colnames(fr.layout) <- paste0("FDG_", 1:2)
  rownames(fr.layout) <- genes.use
  fr.layout$major.modules <- module.DF[rownames(fr.layout), "module"]
  fr.layout$minor.modules <- module.DF[rownames(fr.layout), "minor.module.full"]
  data.nodes <- fr.layout
  data.links <- as.data.frame(A)
  data.links <- data.links %>% mutate(from = rownames(.), .before=1)
  data.links <- data.links %>% pivot_longer(cols = 2:ncol(.), names_to = "to", values_to = "weight")
  data.links <- subset(data.links, weight > 0)
  data.links <- subset(data.links, from != to)
  data.links$x1 <- data.nodes[data.links$from, "FDG_1"]
  data.links$y1 <- data.nodes[data.links$from, "FDG_2"]
  data.links$x2 <- data.nodes[data.links$to, "FDG_1"]
  data.links$y2 <- data.nodes[data.links$to, "FDG_2"]
  if (!"module.DR" %in% names(object@misc[[module.source]])) {
    object@misc[[module.source]]$module.DR <- list("FDG" = list())
  }
  dr.list <- list("nodes" = data.nodes, "links" = data.links)
  if (is.null(on.module)) {
    object@misc[[module.source]]$module.DR$FDG[["all"]] <- dr.list
  } else {
    object@misc[[module.source]]$module.DR$FDG[[on.module]] <- dr.list
  }
  return(object)
}


#' Visualize modular gene interactions using a dimensionality reduction plot
#'
#' This function creates a plot displaying the connections and relationships between genes
#' in major modules, visualized according to a chosen dimensionality reduction technique.
#'
#' @param object A Seurat object containing module dimensionality reduction data.
#' @param reduction A character string specifying which dimensionality reduction technique
#' was used. Default is "FDG" for Force-Directed Graph.
#' @param group.by A character string indicating the field in the `data.nodes`
#' dataframe to use for coloring the points. Default is "major.modules".
#' @param slot A character string indicating the specific module slot to use for the plot.
#' Default is "all".
#' @param module.source A character string specifying the source of the module data within
#' the Seurat object. Defaults to "GeneCooc".
#' @param pt.size A numeric value specifying the size of the points in the plot. Default is 1.
#' @param pt.border A numeric value specifying the size of the point borders. If NULL,
#' no borders will be added to the points. Default is NULL.
#' @param pt.colors A character vector of colors to use for the points. If NULL,
#' a predefined color scale is used. Default is NULL.
#' @param line.color A character string specifying the color of the lines connecting points.
#' Default is "grey".
#' @param line.alpha A numeric value specifying the transparency of the lines connecting
#' points. Default is 0.01.
#' @return A `ggplot` object representing the dimensionality reduction plot.
#' @importFrom ggplot2 ggplot aes geom_segment geom_point scale_linewidth_continuous
#' ggtitle xlab ylab guides theme element_blank element_text guide_legend
#' @importFrom tidydr theme_dr
#' @export
ModuleDimPlot <- function(
    object,
    reduction="FDG",
    group.by="major.modules",
    slot="all",
    module.source="GeneCooc",
    pt.size = 1,
    pt.border=NULL,
    pt.colors=NULL,
    line.color="grey",
    line.alpha=0.01
){
  module.DR.list <- FetchModuleDR(object, reduction = reduction, slot = slot, module.source = module.source)
  data.nodes <- module.DR.list$nodes
  data.links <- module.DR.list$links
  x <- paste0(reduction, "_1")
  y <- paste0(reduction, "_2")
  plot <- ggplot(data.nodes, aes(get(x), get(y))) +
    geom_segment(inherit.aes = F, data = data.links,
                 mapping = aes(x=x1, y=y1, xend=x2, yend=y2, linewidth=weight),
                 alpha=line.alpha, color = line.color, show.legend = F)
  if (!is.null(pt.border)) {
    plot <- plot + geom_point(size = pt.size + pt.border, color = "black")
  }
  plot <- plot +
    geom_point(aes(color = get(group.by)), size = pt.size) +
    scale_linewidth_continuous(breaks = seq(0,1,.2), range = range(0,1))
  if (is.null(pt.colors)){
    plot <- plot + ggsci::scale_color_d3("category20")
  } else {
    plot <- plot + scale_color_manual(values = pt.colors)
  }
  if (slot == "all") {
    title <- "Gene co-occurrance network"
  } else {
    title <- glue::glue("Gene co-occurrance network ({slot})")
  }
  plot <- plot +
    ggtitle(title) + xlab(x) + ylab(y) +
    guides(color = guide_legend(title = "Gene module", override.aes = list(size = 2, alpha=1))) +
    tidydr::theme_dr() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = .5, face = "bold"))
  plot
}
