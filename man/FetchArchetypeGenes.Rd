% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getters_and_setters.R
\name{FetchArchetypeGenes}
\alias{FetchArchetypeGenes}
\title{Fetch Archetype Genes from Seurat Object}
\usage{
FetchArchetypeGenes(object, modules, module.source = "GeneCooc")
}
\arguments{
\item{object}{A Seurat object.}

\item{modules}{A vector of module names specifying which modules to consider when retrieving
archetype genes.}

\item{module.source}{A string specifying the key or named list within the Seurat object from
which to retrieve the module information. The default value is "GeneCooc".}
}
\value{
A named vector of genes.
}
\description{
Retrieves a list of archetype genes from a Seurat object based on the specified module source
and modules. This function returns a vector of gene names that are marked as archetype genes
within the specified modules in the given Seurat object.
}
