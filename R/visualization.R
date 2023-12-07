#' Visualize Gene Modules on UMAP Plot
#'
#' Creates a UMAP plot to display gene modules as colored points. The plot can be customized with
#' different point sizes, transparencies, and an option to highlight specific genes.
#'
#' @param object A Seurat object that contains UMAP coordinates and gene module data.
#' @param pt.size Numeric size of the points in the UMAP plot.
#' @param pt.alpha Numeric value representing the transparency of the points in the UMAP plot.
#' Values range from 0 (transparent) to 1 (opaque).
#' @param text.size Numeric size of the text labels for highlighted genes.
#' @param genes.highlight A vector of gene names to highlight on the UMAP plot. Highlighted
#' genes will be marked with different point shapes and labeled with text.
#' @param colors A named vector of colors to be used for coloring the gene modules. The names
#' should correspond to the unique module names.
#' @param module.source The name of the source from which to fetch the gene module data.
#' Defaults to "GeneCooc".
#'
#' @return A ggplot2 object containing the UMAP plot.
#'
#' @export
#'
ModuleUMAPPlot <- function(
    object,
    pt.size = 1,
    pt.alpha = 1,
    text.size = 3,
    genes.highlight = NULL,
    colors = NULL,
    module.source = "GeneCooc"
){
  mods <- FetchModuleDF(object, module.source = module.source)
  mods <- subset(mods, !is.na(UMAP_1))
  p <- ggplot(mods, aes(UMAP_1, UMAP_2, color = module) ) +
    geom_point(size = pt.size, alpha = pt.alpha) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  if (!is.null(genes.highlight)) {
    p <- p +
      geom_point(
        inherit.aes = F,
        data = subset(mods, gene.name %in% genes.highlight),
        aes(UMAP_1, UMAP_2),
        size = pt.size+0.5,
        shape = 1,
        color = "black") +
      ggrepel::geom_text_repel(
        inherit.aes = F,
        data = subset(mods, gene.name %in% genes.highlight),
        aes(UMAP_1, UMAP_2, label = gene.name),
        color="blue",
        size = text.size)
  }
  p + theme_bw(base_size = 15) + theme(legend.title = element_blank())
}


#' Visualize Specific Gene Sub-Modules on UMAP Plot
#'
#' Generates a UMAP plot that highlights specific gene sub-modules within a Seurat object.
#' Additional genes may be highlighted with distinct markers and labels for further emphasis.
#'
#' @param object A Seurat object that contains UMAP coordinates and gene module data.
#' @param pt.size Numeric size of the points in the UMAP plot.
#' @param pt.alpha Numeric value representing the transparency of the points in the UMAP plot.
#' Values range from 0 (transparent) to 1 (opaque).
#' @param text.size Numeric size of the labels for highlighted genes on the UMAP plot.
#' @param module.use A character vector indicating which sub-modules should be highlighted on
#' the UMAP plot. This parameter is required.
#' @param genes.highlight A vector of gene names to highlight on the plot.
#' @param colors A named vector of colors to be used for the sub-modules. The names should
#' correspond to the unique names of the sub-modules you wish to color.
#' @param module.source The name of the source from which to fetch the gene module data.
#' Defaults to "GeneCooc".
#'
#' @return A `ggplot2` object containing the UMAP plot with distinct visual emphasis on
#' the selected gene sub-modules.
#'
#' @export
#'
SubModuleUMAPPlot <- function(
    object,
    pt.size = 1,
    pt.alpha = 1,
    text.size = 3,
    module.use = NULL,
    genes.highlight = NULL,
    colors = NULL,
    module.source = "GeneCooc"
){
  if (is.null(module.use)) {
    stop("Parameter `module.use` were not provided!")
  }
  mods <- FetchModuleDF(object, module.source = module.source)
  mods <- subset(mods, !is.na(UMAP_1))
  mods.filter <- subset(mods, module %in% module.use)
  mods.bg <- subset(mods, !module %in% module.use)
  p <- ggplot(mods.filter, aes(UMAP_1, UMAP_2, color = minor.module.full) ) +
    geom_point(inherit.aes = F, data = mods.bg, aes(UMAP_1, UMAP_2),
               size = pt.size, alpha = pt.alpha, color = 'grey') +
    geom_point(size = pt.size, alpha = pt.alpha) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))
  if(!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  if (!is.null(genes.highlight)) {
    p <- p +
      geom_point(
        inherit.aes = F,
        data = subset(mods.filter, gene.name %in% genes.highlight),
        aes(UMAP_1, UMAP_2),
        size = pt.size+0.5,
        shape = 1,
        color = "black") +
      ggrepel::geom_text_repel(
        inherit.aes = F,
        data = subset(mods.filter, gene.name %in% genes.highlight),
        aes(UMAP_1, UMAP_2, label = gene.name, color = minor.module.full),
        color="blue",
        size = text.size)
  }
  p + theme_bw(base_size = 15) + theme(legend.title = element_blank())
}

