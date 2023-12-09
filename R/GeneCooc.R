#' @include getters_and_setters.R
NULL

#' Calculate Gene Rankings
#'
#' This function calculates gene rankings using Multiple Correspondence Analysis (MCA) implemented in CelliD.
#'
#' @param object A Seurat object.
#' @param min.expr.cells The minimum number of cells in which a gene must be expressed
#' to be considered. Defaults to 20. If min.expr.pct is specified, min.expr.cells is
#' adjusted to be the minimal between min.expr.cells and ncol(object)*min.expr.pct.
#' @param min.expr.pct The minimum percentage of cells in which a gene must be expressed
#' to be considered. Defaults to 0.001.
#' @param nfeatures The number of features to choose when finding variable features.
#' Defaults to 8000.
#' @param features An optional vector of specific features to use for the analysis.
#' If NULL, features will be automatically selected based on expression criteria.
#' Defaults to NULL.
#' @param ndim.mca The number of dimensions to consider in the MCA. Defaults to 30.
#' @param use.variable.features Logical, indicating whether to only consider variable
#' features. If TRUE, the function searches for variable features using the `FindVariableFeatures`
#' function from Seurat. Defaults to TRUE.
#' @param batch.name An optional character string representing the batch variable
#' to account for batch effects when selecting variable features. If NULL, batch effects are not considered during
#' the variable features selection. Defaults to NULL.
#' @param module.source A character string indicating where to save results of `GeneCooc`. Default is "GeneCooc".
#' @param assay A character string specifying which assay to use from the Seurat
#' object. Defaults is "RNA".
#' @return A modified Seurat object inclusive of the `features.use` and `gene.rankings` in its
#' `misc` slot in the list named by `module.source`.
#' @import Seurat
#' @export
CalGeneRankings <- function(
    object,
    min.expr.cells=20,
    min.expr.pct=0.001,
    nfeatures=8000,
    features=NULL,
    ndim.mca=30,
    use.variable.features=TRUE,
    batch.name=NULL,
    module.source="GeneCooc",
    assay="RNA"
){
  ## get features
  expr.mat <- GetAssayData(object, assay = 'RNA', slot = 'counts')
  if (is.null(features)) {
    min.expr.cells <- min(min.expr.cells, ncol(expr.mat) * min.expr.pct)
    expr.in.cells <- Matrix::rowSums(expr.mat > 0)
    features.use <- names(expr.in.cells)[expr.in.cells >= min.expr.cells]
  } else {
    features.use <- intersect(rownames(expr.mat), features)
  }
  if (length(features.use) == 0) {
    stop('No features error.')
  }
  ## find variable genes
  if(use.variable.features) {
    if (is.null(batch.name)) {
      object <- FindVariableFeatures(object, nfeatures = nfeatures, selection.method = "vst")
    } else {
      object.list <- SplitObject(object, split.by = batch.name)
      for (i in seq_along(object.list)) {
        object.list[[i]] <- FindVariableFeatures(object, nfeatures = nfeatures, selection.method = "vst")
      }
      VariableFeatures(object) <- SelectIntegrationFeatures(object.list = object.list, nfeatures = nfeatures)
    }
    features.use <- intersect(features.use, VariableFeatures(object))
  }
  ## run MCA
  object <- CelliD::RunMCA(object, features = features.use, nmcs = ndim.mca)
  gene.rankings <- CelliD::GetCellGeneRanking(object, dims = 1:ndim.mca)
  ## set results to object@misc slot
  object@misc[[module.source]] <- list(
    "features.use" = features.use,  # vector (character)
    "gene.rankings" = gene.rankings # list
  )
  return(object)
}



#' Calculate Affinity Matrix
#'
#' This function calculates the affinity (co-occurrence) matrix for a given set of features
#' within the provided Seurat object, based on the top K gene rankings. The affinity matrix represents
#' the co-occurrence of genes across gene rankings. Genes that occur together more frequently
#' receive higher scores in the affinity matrix. Genes with a frequency below a minimum threshold
#' are trimmed from the analysis.
#'
#' @param object A Seurat object containing `gene.rankings` and `features.use` information within the
#' `Misc` slot.
#' @param K An integer value specifying the number of top genes to consider for each ranking.
#' Default is 500.
#' @param min.freq An integer value indicating the minimum frequency a gene must have to be
#' included in the final affinity matrix. Genes appearing less frequently are dropped. Default
#' is 10.
#' @param module.source A character string indicating where to load tmp results and save final results
#' of `GeneCooc`. Default is "GeneCooc".
#'
#' @return A modified Seurat object inclusive of the calculated affinity matrix in its
#' `Misc` slot in the list named by `module.source`.
#'
#' @import Seurat
#' @export
CalAffinityMatrix <- function(object, K=500, min.freq=10, module.source="GeneCooc"){
  ## get top-K genes
  message("Fetching top K genes ...")
  gene.rankings <- Misc(object)[[module.source]]$gene.rankings
  features.use <- Misc(object)[[module.source]]$features.use
  sentences <- lapply(gene.rankings, function(xx) names(head(xx, K)))
  ## init affinity matrix
  A <- matrix(0 , nrow = length(features.use), ncol = length(features.use))
  rownames(A) <- features.use
  colnames(A) <- rownames(A)
  ## calculate affinity (co-occurrence) matrix
  message("Calculating affinity matrix ...")
  batch.size <- round(length(sentences) / 10)
  process <- 0
  for(i in seq_along(sentences)) {
    if (i %% batch.size == 0) {
      process <- process + 1
      message(glue::glue("Processing {process}/10 ..."))
    }
    s <- sentences[[i]]
    A[s, s] <- A[s, s] + 1
  }
  ## trim genes
  message("Trimming genes ...")
  idx <- diag(A) >= min.freq
  A <- A[idx, idx]
  n.dropped.genes <- table(idx)['FALSE']
  message(glue::glue("{n.dropped.genes} of {nrow(A)} genes appeared less than {min.freq} times were dropped."))
  message("Normalizing ...")
  A <- apply(A, 2, function(xx) xx / max(xx))
  A <- (A + t(A)) / 2
  ## write results into Seurat object
  object@misc[[module.source]]$affinity.matrix <- A
  return(object)
}


#' Find Gene Modules
#'
#' This function identifies gene modules based on an affinity matrix within a given Seurat object.
#' It uses the Louvain algorithm for primary clustering and a dynamic tree cut for sub-module detection.
#'
#' @param object A Seurat object containing gene co-occurrence data in an affinity matrix within
#' its `misc` slot.
#' @param k Integer. Specifies the number of nearest neighbors to be used in the construction of
#' the SNN (shared nearest neighbor) graph.
#' @param resolution A double. Resolution parameter for the Louvain clustering algorithm. Controls
#' the granularity of the clustering.
#' @param min.module.size Integer. The minimum size of the modules to be detected in the dynamic
#' tree cut.
#' @param module.source A character string indicating where to load tmp results and save final results
#' of `GeneCooc`. Default is "GeneCooc".
#'
#' @return Seurat object with the gene module information computed. The results will be written
#' into the object's `misc` slot under `GeneCooc`.
#'
#' @export
FindModules <- function(object, k=50, resolution=0.1, min.module.size=10, module.source="GeneCooc"){
  ## fetch data
  A <- Misc(object)[[module.source]]$affinity.matrix
  dissM <- 1 - A
  ## find major modules
  g <- Seurat::FindNeighbors(as.dist(dissM), k.param = 50)
  g <- igraph::graph_from_adjacency_matrix(adjmatrix = g$snn,
                                           mode = "undirected",
                                           weighted = TRUE)
  ## Louvain cluster
  clusters <- igraph::cluster_louvain(g, resolution = resolution)
  mods <- data.frame(
    row.names = rownames(A),
    gene.name = rownames(A),
    module = paste0(glue::glue("{module.source}-M"), clusters$membership)
  )

  ## find sub-modules in each major module: dynamic tree cut
  mods$minor.module <- NA
  module.names <- sort(unique(mods$module))
  for (mn in module.names) {
    genes <- subset(mods, module == mn)$gene.name
    D <- dissM[genes, genes]
    gene.tree <- stats::hclust(d = as.dist(D), method = 'average')
    dynamic.mods = dynamicTreeCut::cutreeDynamic(dendro = gene.tree, distM = D,
                                                 minClusterSize = min.module.size)
    mods[gene.tree$labels, 'minor.module'] <- dynamic.mods
  }
  mods$minor.module.full <- paste0(mods$module, '-', mods$minor.module)
  ## write results into Seurat object
  object@misc[[module.source]]$gene.module <- mods
  return(object)
}


#' Trim Gene Modules Based on Archetype Analysis
#'
#' Performs trimming of gene sub-modules within the given object using archetype analysis. The
#' function evaluates each sub-module based on archetypes and trims them according to specified
#' cutoff values.
#'
#' @param object A Seurat object containing gene module data and an affinity matrix within its
#' `misc` slot.
#' @param archetype.score.cutoff A numeric threshold for the archetype score. Genes with an
#' archetype score below this cutoff will not be considered archetypes and may be trimmed
#' accordingly.
#' @param delta.cutoff A numeric threshold for the difference in archetype scores. If the difference
#' between the scores of the most expressive genes (archetypes) is less than this value, they may be
#' excluded from the module.
#' @param module.source A character string indicating where to load tmp results and save final results
#' of `GeneCooc`. Default is "GeneCooc".
#'
#' @return Seurat object with updated sub-module information, including flags to indicate which
#' genes are kept and which ones are identified as archetypes for their respective sub-modules.
#'
#' @importFrom magrittr `%>%`
#'
#' @export
TrimModules <- function(object, archetype.score.cutoff=0.5, delta.cutoff=0.2, module.source="GeneCooc") {
  ## fetch data
  mods <- Misc(object)[[module.source]]$gene.module
  A <- Misc(object)[[module.source]]$affinity.matrix
  ## trim sub-modules via archetypes analysis
  module.names <- sort(unique(mods$minor.module.full))
  mods$is.kept <- FALSE
  mods$is.archetype <- FALSE
  for (module.name in module.names[!endsWith(module.names, '-0')]) {
    message(glue::glue('processing {module.name} ...'))
    genes <- subset(mods, minor.module.full == module.name)$gene.name
    a <- A[genes, genes]
    ## archetype analysis
    ## while loop to aviod invalid results from archetypes()
    archetype.mat <- NULL
    while (!is.matrix(archetype.mat)) {
      archetype.mat <- archetypes::archetypes(a, k = 2, verbose = FALSE)$archetypes
    }
    ## trim sub-modules
    delta <- abs(apply(archetype.mat, 2, diff))
    gene.ats <- apply(archetype.mat, 2, which.max)
    max.ind <- names(table(gene.ats))[which.max(table(gene.ats))] %>% as.integer()
    gene.ats[archetype.mat[max.ind, ] < archetype.score.cutoff] <- 0
    gene.ats[delta < delta.cutoff] <- 0
    genes.kept <- genes[gene.ats == max.ind]
    if (length(genes.kept) > 0) {
      mods[genes.kept, ]$is.kept <- TRUE
    }
    ## archetype gene for each sub-module
    mods[names(which.max(archetype.mat[max.ind, ])), ]$is.archetype <- TRUE
  }
  ## write results into Seurat object
  object@misc[[module.source]]$gene.module <- mods
  return(object)
}


#' Run UMAP on Gene Modules
#'
#' This function runs UMAP (Uniform Manifold Approximation and Projection) on gene modules within
#' a given Seurat object. It supports both supervised and unsupervised modes.
#'
#' @param object A Seurat object containing gene module data and an affinity matrix within its
#' `misc` slot.
#' @param exclude.trimmed Logical, indicating whether to exclude trimmed genes from UMAP analysis.
#' If `TRUE`, genes that are not marked as kept will be omitted.
#' @param supervised Logical, controls the mode of UMAP execution. If `FALSE`, the UMAP runs in
#' unsupervised mode. If `TRUE`, it utilizes gene module labels to steer the embedding in a
#' supervised fashion.
#' @param module.source A character string indicating where to load tmp results and save final results
#' of `GeneCooc`. Default is "GeneCooc".
#' @param ... Additional arguments passed on to the uwot::umap function.
#'
#' @return Seurat object with updated gene module information including UMAP coordinates.
#'
#' @export
#'
RunModuleUMAP <- function(object, exclude.trimmed=TRUE, supervised=FALSE, module.source="GeneCooc", ...) {
  ## fetch data
  mods.origin <- Misc(object)[[module.source]]$gene.module
  mods <- Misc(object)[[module.source]]$gene.module
  A <- Misc(object)[[module.source]]$affinity.matrix
  if (exclude.trimmed){
    mods <- subset(mods, is.kept)
    A <- A[rownames(mods), rownames(mods)]
  }
  dissM <- 1 - A
  diag(dissM) <- 0
  ## unsupervised UMAP
  if (supervised) {
    embeddings <-  uwot::umap(X = as.dist(dissM), y = mods$minor.module.full, ...)
  } else {
    embeddings <-  uwot::umap(X = as.dist(dissM), ...)
  }
  mods.origin$UMAP_1 <- NA
  mods.origin$UMAP_2 <- NA
  mods.origin[rownames(embeddings), "UMAP_1"] <- embeddings[, 1]
  mods.origin[rownames(embeddings), "UMAP_2"] <- embeddings[, 2]
  ## write results into Seurat object
  object@misc[[module.source]]$gene.module <- mods.origin
  return(object)
}


#' Calculate Module Scores for Cells in a Seurat Object
#'
#' Computes scores for gene modules in each cell based on multiple correspondence analysis (MCA)
#' embedding and the association of genes to modules. It first filters out genes that have been
#' trimmed, then calculates cell-to-gene distances, followed by aggregating these distances to
#' compute module scores for each cell.
#'
#' @param object A Seurat object containing gene module data in its `misc` slot,
#' and embeddings from MCA.
#' @param ndim.mca The number of dimensions to use from the multiple correspondence analysis
#' (MCA) embeddings for calculating the module scores. Default is 30 dimensions.
#' @param min.size Minimal size of the gene set for calculating the module score. Default: 10.
#' @param module.source A character string indicating where to load tmp results and save final results
#' of `GeneCooc`. Default is "GeneCooc".
#'
#' @return Seurat object updated with a new assay slot named `module.score` containing the
#' calculated scores for each module and cell.
#'
# TODO calculate module scores for a given gene module list
#' @export
#'
CalModuleScore <- function(object, ndim.mca=30, min.size=10, module.source="GeneCooc") {
  ## fetch data
  mods <- Misc(object)[[module.source]]$gene.module
  mods <- subset(mods, is.kept) ## drop the trimmed genes
  X <- Embeddings(object, reduction = "mca")[, 1:ndim.mca]
  Y <- Loadings(object, reduction = "mca")[rownames(mods), 1:ndim.mca]
  ## calculate gene to cell distance
  Z <- proxy::dist(X, Y, method = "Euclidean")
  Z <- scale(t(1/Z))
  rownames(Z) <- rownames(Y) # rows: genes
  colnames(Z) <- rownames(X) # cols: cells
  ## calculate module score
  module.list <- FetchModuleList(object, module.source = module.source, module.type = "both")
  module.size <- sapply(module.list, length)
  module.used <- names(module.size)[module.size >= min.size]
  module.list <- module.list[module.used]
  M <- sapply(module.list, function(xx) {
    as.numeric(rownames(Z) %in% xx)
  })
  ## module score
  ## Z(n.genes, n.cells), M(n.genes, n.modules), score(n.cells, n.modules)
  score <- crossprod(Z, M) ## t(Z) %*% M
  score <- scale(score)
  object[[module.source]] <- CreateAssayObject(data = t(score))
  return(object)
}

