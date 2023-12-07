FetchModuleDF <- function(object, module.source="GeneCooc"){
  return(Misc(object)[[module.source]]$gene.module)
}


FetchAffinityMatrix <- function(object, module.source="GeneCooc"){
  return(Misc(object)[[module.source]]$affinity.matrix)
}


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


