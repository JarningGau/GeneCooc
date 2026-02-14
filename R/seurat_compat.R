#' @keywords internal
.GeneCoocSeuratMajorVersion <- function() {
  as.integer(utils::packageVersion("Seurat")[1])
}

#' @keywords internal
.GeneCoocGetCounts <- function(object, assay = "RNA") {
  if (.GeneCoocSeuratMajorVersion() >= 5) {
    return(Seurat::GetAssayData(object, assay = assay, layer = "counts"))
  }
  Seurat::GetAssayData(object, assay = assay, slot = "counts")
}

#' @keywords internal
.GeneCoocFindNeighbors <- function(dissM, k = 50) {
  g <- Seurat::FindNeighbors(as.dist(dissM), k.param = k)
  if ("snn" %in% names(g)) {
    return(g$snn)
  }
  snn.name <- grep("_snn$", names(g), value = TRUE)
  if (length(snn.name) > 0) {
    return(g[[snn.name[[1]]]])
  }
  stop("Could not find the SNN graph in Seurat::FindNeighbors output.")
}
