#' @include getters_and_setters.R
NULL

.ModulePairDist <- function(A, module.list, module.1, module.2) {
  genes.1 <- module.list[[module.1]]
  genes.2 <- module.list[[module.2]]
  d <- (mean(A[genes.1, genes.1]) + mean(A[genes.2, genes.2])) / 2 - mean(A[genes.1, genes.2])
  return(d)
}


ModuleDist <- function(object, module.source="GeneCooc", major.module=NULL) {
  A <- FetchAffinityMatrix(object, module.source)
  module.list <- FetchModuleList(object, module.source, module.type = "minor")
  mods <- FetchModuleDF(object, module.source)
  mods <- subset(mods, is.kept) ## drop the trimmed genes
  if (!is.null(major.module)) {
    mods <- subset(mods, module == major.module)
  }
  minor.modules <- sort(unique(mods$minor.module.full))
  n.modules <- length(minor.modules)
  D <- matrix(1, nrow = n.modules, ncol = n.modules)
  rownames(D) <- minor.modules
  colnames(D) <- minor.modules
  module.pairs <- combn(minor.modules, m = 2, simplify = F)
  module.pairs <- lapply(module.pairs, sort)
  for(m in module.pairs) {
    d <- .ModulePairDist(A, module.list, m[1], m[2])
    D[m[1], m[2]] <- d
    D[m[2], m[1]] <- d
  }
  diag(D) <- 0
  return(D)
}


NextRelatedModulePairs <- function(object, module.source="GeneCooc", major.module=NULL, do.plot=TRUE){
  paired.dist <- ModuleDist(object, module.source, major.module)
  if (nrow(paired.dist) == 2) {
    paired.module <- rownames(paired.dist)
  } else {
    hc <- hclust(as.dist(paired.dist), method = "average")
    if (do.plot) {
      plot(hc)
    }
    idx <- hc$merge[1, ]
    paired.module <- hc$labels[-idx]
  }
  return(paired.module)
}


EstimateAccuracy <- function(object, module.1, module.2, module.source="GeneCooc", do.plot=TRUE) {
  A <- FetchAffinityMatrix(object, module.source)
  mods <- FetchModuleDF(object, module.source)
  module.list <- FetchModuleList(object, module.source, module.type = "minor")
  genes.1 <- module.list[[module.1]]
  genes.2 <- module.list[[module.2]]
  all.genes <- c(genes.1, genes.2)
  a <- A[all.genes, all.genes]
  archetype.mat <- NULL
  while (!is.matrix(archetype.mat)) {
    archetype.mat <- archetypes::archetypes(a, k = 2, verbose = FALSE)$archetypes
  }
  rownames(archetype.mat) <- paste("Archetype", 1:2)
  gene.labels <- ifelse(all.genes %in% genes.1, 1, 2)
  top.annotation <- mods[all.genes, "minor.module.full", drop=F]
  if (do.plot) {
    pheatmap::pheatmap(archetype.mat, annotation_col = top.annotation, show_colnames = F,
                       show_rownames = T, cluster_rows = F, cutree_cols = 2,
                       clustering_distance_cols = "euclidean",
                       clustering_method = "average")
  }
  hc <- hclust(dist(t(archetype.mat)), method = "average")
  clusters <- dendextend::cutree(hc, k=2)
  confusion.matrix <- as.matrix(table(clusters, gene.labels))
  confusion.matrix <- confusion.matrix[order(confusion.matrix[1, ], decreasing = T), ]
  sum(diag(confusion.matrix)) / sum(confusion.matrix)
}


MergeModules <- function(object, module.1, module.2, module.source="GeneCooc") {
  message(glue::glue("Merging {module.2} to {module.1}"))
  mods <- FetchModuleDF(object, module.source)
  genes.1 <- rownames(subset(mods, minor.module.full == module.1))
  genes.2 <- rownames(subset(mods, minor.module.full == module.2))
  if (!"minor.module.full.before.merge" %in% colnames(mods)) {
    mods$minor.module.full.before.merge <- mods$minor.module.full
  }
  mods[genes.2, "minor.module.full"] <- module.1
  object@misc[[module.source]]$gene.module <- mods
  return(object)
}

.StatMinorModules <- function(object, module.source="GeneCooc", major.module=NULL) {
  mods <- FetchModuleDF(object, module.source)
  mods <- subset(mods, is.kept) ## drop the trimmed genes
  if (!is.null(major.module)) {
    mods <- subset(mods, module == major.module)
  }
  n.minor.modules <- length(unique(mods$minor.module.full))
  return(n.minor.modules)
}

#' Auto Merge Modules Function
#'
#' @description
#' This function iteratively merges minor modules within major modules based on an accuracy
#' threshold until all minor modules are processed or accuracy is no longer improved beyond the
#' threshold.
#'
#' @param object A Seurat object.
#' @param acc.threshold The accuracy threshold for merging modules, two modules will merge
#' if the classification accuracy less than this value; defaults to 0.9.
#' @param min.size Minimal size of gene module for module merging. If one of two gene modules
#' to be merged less than this value, it will merge them without consideration of `acc.threshold` Default: 2.
#' @param module.source A character string indicating where to save results of `GeneCooc`. Default is "GeneCooc".
#'
#' @return A modified Seurat object.
#'
#' @export
#'
AutoMergeModules <- function(object, acc.threshold=0.9, min.size=2, module.source="GeneCooc"){
  ## fetch data
  module.list <- FetchModuleList(object, module.source, module.type = "major")
  for (major.module in names(module.list)) {
    message("===================")
    message(glue::glue("Processing {major.module} ..."))
    n.minor.modules <- .StatMinorModules(object, module.source, major.module)
    if (n.minor.modules > 1) {
      acc <- 0
      while(TRUE) {
        paired.module <- NextRelatedModulePairs(object, module.source, major.module, do.plot = F)
        acc <- EstimateAccuracy(object, paired.module[1], paired.module[2], module.source, do.plot = F)
        message(glue::glue("Accuracy of {paired.module[1]} and {paired.module[2]}: {round(acc,3)}"))
        if (acc < acc.threshold || length(paired.module[1]) <= min.size || length(paired.module[2]) <= min.size) {
          object <- MergeModules(object, paired.module[1], paired.module[2], module.source)
        }
        n.minor.modules <- .StatMinorModules(object, module.source, major.module)
        message(glue::glue("{n.minor.modules} minor modules left in {major.module}"))
        if (acc >= acc.threshold || n.minor.modules < 2) {
          message(glue::glue("Finished: {major.module}"))
          message("===================")
          break
        }
      }
    }
  }
  return(object)
}
