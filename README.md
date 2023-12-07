
# GeneCooc

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![](https://img.shields.io/badge/devel%20version-0.0.0.9000-green.svg)](https://github.com/jarninggau/GeneCooc)

GeneCooc is an ultra-fast R package for gene co-expression module
discovery facilitates cell state identification.

## Installation

You can install the development version of `GeneCooc` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# install CelliD
devtools::install_github("RausellLab/CelliD")
# install GeneCooc
devtools::install_github("JarningGau/GeneCooc")
```

Note: `GeneCooc` is compitable with Seurat V4, we did not test the codes
on Seurat V5.

## Quick start

``` r
## 1. Calculate the gene to cell distance and rank the genes by this distance for each cell.
seu <- CalGeneRankings(seu, min.expr.cells = 100)
## 2. Calculate the gene affinity matrix. The gene affinity is measured by cooccurance ratio.
seu <- CalAffinityMatrix(seu, K = 200, min.freq = 10)
## 3. Find modules. The major modules are divided by louvain cluster on gene-gene coexpression 
##    graph defined by gene affinity matrix. Then minor modules are divided using the dynamic 
##    tree cut on a hierarchical tree for each major module.
seu <- FindModules(seu)
## 4. Trim the minor modules by archytype analysis.
seu <- TrimModules(seu)
## 5. Merge the similary minor modules automaticall。
seu <- AutoMergeModules(seu, acc.threshold = 0.99)
## 6. Scoring each gene module.
seu <- CalModuleScore(seu)
```
