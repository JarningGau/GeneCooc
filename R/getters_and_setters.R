#' Fetch Module Data Frame From Seurat Object
#'
#' @description
#' Retrieves a data frame of gene modules from a Seurat object based on the specified
#' `module.source`. This is typically used to access a specific assay or custom data that has been
#' added to the Seurat object under a particular key.
#'
#' @param object A Seurat object.
#' @param module.source A character string indicating where to fetch data. Default is "GeneCooc".
#' @param module.size Only returns the modules with sizes more than `module.size`. Default: 0.
#' @param is.trimmed Only return trimmed modules. Default: FALSE.
#'
#' @return A data frame of gene modules from the specified source within the object.
#'
#' @export
FetchModuleDF <- function(object, module.source="GeneCooc", module.size=0, is.trimmed=FALSE){
  mods <- Misc(object)[[module.source]]$gene.module
  if (is.trimmed && "is.kept" %in% colnames(mods)) {
    mods <- subset(mods, is.kept)
  }
  minor.modules <- sort(unique(mods$minor.module.full))
  minor.module.list <- lapply(minor.modules, function(xx) {
    subset(mods, minor.module.full == xx)$gene.name
  })
  names(minor.module.list) <- minor.modules
  module.sizes <- sapply(minor.module.list, length)
  minor.module.list <- minor.module.list[module.sizes >= module.size]
  mods <- subset(mods, minor.module.full %in% names(minor.module.list))
  return(mods)
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
#' @param is.trimmed Only return trimmed modules. Default: TRUE
#' @param module.size Only returns the modules with sizes more than `module.size`. Default: 0.
#' @param module.type A string specifying the type of gene modules to include in the output list:
#' "both" for both major and minor modules, "major" for major modules only, or "minor" for minor modules
#' only. The default value is "both".
#' @param module.source A string specifying the key or named list within the Seurat object from
#' which to retrieve the module information. The default value is "GeneCooc".
#'
#' @return A list of gene names associated with the specified module source and type in the Seurat
#' object. The list is organized based on the major and minor module categories, if applicable.
# TODO: unbound names will cause bugs, solve these unbound names.
#' @export
FetchModuleList <- function(object, is.trimmed=TRUE, module.size=0, module.type="both", module.source="GeneCooc") {
  mods <- FetchModuleDF(object, module.source = module.source, module.size = module.size, is.trimmed = is.trimmed)
  # if ("is.kept" %in% colnames(mods)) {
  #   mods <- subset(mods, is.kept) ## drop the trimmed genes
  # }
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
  module.sizes <- sapply(module.list, length)
  module.list <- module.list[module.sizes >= module.size]
  return(module.list)
}


#' Fetch Archetype Genes from Seurat Object
#'
#' @description
#' Retrieves a list of archetype genes from a Seurat object based on the specified module source
#' and modules. This function returns a vector of gene names that are marked as archetype genes
#' within the specified modules in the given Seurat object.
#'
#' @param object A Seurat object.
#' @param modules A vector of module names specifying which modules to consider when retrieving
#' archetype genes.
#' @param module.source A string specifying the key or named list within the Seurat object from
#' which to retrieve the module information. The default value is "GeneCooc".
#'
#' @return A named vector of genes.
#'
#' @export
FetchArchetypeGenes <- function(object, modules, module.source="GeneCooc") {
  mods <- FetchModuleDF(object, module.source = module.source)
  mods <- subset(mods, is.archetype)
  mods <- subset(mods, minor.module.full.before.merge %in% modules)
  genes <- mods$gene.name
  names(genes) <- modules
  return(genes)
}


#' Rescue Specified Modules
#'
#' @description
#' This function "rescues" certain gene modules in a Seurat object by marking them as kept,
#' effectively including them for further analysis or consideration. The specified modules
#' are altered to have their `is.kept` status set to `TRUE` within the module source data.
#'
#' @param object A Seurat object containing various types of data, including gene module information.
#' @param modules A character vector of minor module names (`minor.module.full`) to be marked as kept
#'  within the Seurat
#' object.
#' @param module.source A string specifying the name of the list key within the Seurat object's
#' `@misc` slot from which gene module information should be fetched and updated.
#' The default value is "GeneCooc".
#'
#' @return The modified Seurat object with the specified subset of modules marked as "kept".
#'
#' @export
RescueModules <- function(object, modules, module.source="GeneCooc") {
  mods <- FetchModuleDF(object, module.source = module.source)
  mods[mods$minor.module.full %in% modules, "is.kept"] = TRUE
  object@misc[[module.source]]$gene.module <- mods
  return(object)
}


#' Fetch Module Dimensionality Reduction Data
#'
#' This function retrieves dimensionality reduction data for gene modules
#' from a specified slot within a given Seurat object.
#'
#' @param object A Seurat object containing module dimensionality reduction information.
#' @param reduction A character string indicating which type of dimensionality reduction data
#' to fetch. Default is "FDG" for Force-Directed Graph data.
#' @param slot A character string specifying the module slot to retrieve data from.
#' Default is "all" which would typically denote the combined results for all modules.
#' @param module.source A character string specifying the key under which the module
#' dimensionality reduction data is stored within the Seurat object's `misc` slot. Default is "GeneCooc".
#'
#' @return Returns a data structure containing the dimensionality reduction data.
#' The exact structure will depend on how the data is stored within the provided Seurat object.
#' @export
FetchModuleDR <- function(
    object,
    reduction="FDG",
    slot="all",
    module.source="GeneCooc"
){
  return(object@misc[[module.source]]$module.DR[[reduction]][[slot]])
}
