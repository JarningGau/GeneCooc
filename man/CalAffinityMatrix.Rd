% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneCooc.R
\name{CalAffinityMatrix}
\alias{CalAffinityMatrix}
\title{Calculate Affinity Matrix}
\usage{
CalAffinityMatrix(object, K = 500, min.freq = 10, module.source = "GeneCooc")
}
\arguments{
\item{object}{A Seurat object containing \code{gene.rankings} and \code{features.use} information within the
\code{Misc} slot.}

\item{K}{An integer value specifying the number of top genes to consider for each ranking.
Default is 500.}

\item{min.freq}{An integer value indicating the minimum frequency a gene must have to be
included in the final affinity matrix. Genes appearing less frequently are dropped. Default
is 10.}

\item{module.source}{A character string indicating where to load tmp results and save final results
of \code{GeneCooc}. Default is "GeneCooc".}
}
\value{
A modified Seurat object inclusive of the calculated affinity matrix in its
\code{Misc} slot in the list named by \code{module.source}.
}
\description{
This function calculates the affinity (co-occurrence) matrix for a given set of features
within the provided Seurat object, based on the top K gene rankings. The affinity matrix represents
the co-occurrence of genes across gene rankings. Genes that occur together more frequently
receive higher scores in the affinity matrix. Genes with a frequency below a minimum threshold
are trimmed from the analysis.
}
