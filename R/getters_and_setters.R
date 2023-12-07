#' Fetch Module Data Frame From Seurat Object
#'
#' @description
#' Retrieves a data frame of gene modules from a Seurat object based on the specified
#' `module.source`. This is typically used to access a specific assay or custom data that has been
#' added to the Seurat object under a particular key.
#'
#' @param object A Seurat object.
#' @param module.source A character string indicating where to fetch data. Default is "GeneCooc".
#'
#' @return A data frame of gene modules from the specified source within the object.
#'
#' @export
FetchModuleDF <- function(object, module.source="GeneCooc"){
  return(Misc(object)[[module.source]]$gene.module)
}


#' Fetch Affinity Matrix From Seurat Object
#'
#' @description
#' Retrieves an affinity matrix from a Seurat object based on the specified module source.
#'
#' @param object A Seurat object.
#' @param module.source A character string indicating where to fetch data. Default is "GeneCooc".
#'
#' @return An affinity matrix associated with the specified key in the Seurat object's
#' miscellaneous data slot or another designated slot where such data is stored.
#'
#' @export
FetchAffinityMatrix <- function(object, module.source="GeneCooc"){
  return(Misc(object)[[module.source]]$affinity.matrix)
}


#' Fetch Gene Module List From Seurat Object
#'
#' @description
#' Retrieves a gene module list from a Seurat object based on the specified module source and type.
#' This function returns a list of gene names associated with different modules within the provided
#' Seurat object.
#'
#' @param object A Seurat object.
#' @param module.source A string specifying the key or named list within the Seurat object from
#' which to retrieve the module information. The default value is "GeneCooc".
#' @param module.type A string specifying the type of gene modules to include in the output list:
#' "both" for both major and minor modules, "major" for major modules only, or "minor" for minor modules
#' only. The default value is "both".
#'
#' @return A list of gene names associated with the specified module source and type in the Seurat
#' object. The list is organized based on the major and minor module categories, if applicable.
#'
#' @export
FetchModuleList <- function(object, module.source="GeneCooc", module.type="both") {
  mods <- FetchModuleDF(object, module.source)
  mods <- subset(mods, is.kept) ## drop the trimmed genes
  major.modules <- sort(unique(mods$module))
  minor.modules <- sort(unique(mods$minor.module.full))
  major.module.list <- lapply(major.modules, function(xx) {
    subset(mods, module == xx)$gene.name
  })
  minor.module.list <- lapply(minor.modules, function(xx) {
    subset(mods, minor.module.full == xx)$gene.name
  })
  if (module.type == "both") {
    module.list <- c(major.module.list, minor.module.list)
    names(module.list) <- c(major.modules, minor.modules)
  } else if(module.type == "major") {
    module.list <- major.module.list
    names(module.list) <- major.modules
  } else if(module.type == "minor") {
    module.list <- minor.module.list
    names(module.list) <- minor.modules
  } else {
    stop("Invalid parameter: module.type should be one of 'both', 'major' or 'minor'.")
  }
  return(module.list)
}


