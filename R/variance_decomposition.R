#' @importFrom magrittr "%>%" "%<>%"
NULL

#' Perform Variance Decomposition
#'
#' This function performs variance decomposition for a set of variables of interest (`vd.vars`)
#' across genes within the provided Seurat object. Variance decomposition determines the
#' contribution of each variable of interest to the total variance in gene expression. The
#' function uses a linear mixed-effects model to estimate the variance components. The result is
#' a table containing the variance contribution for each variable, as well as the residual
#' variance.
#'
#' @param object A Seurat object containing gene expression data and meta data.
#' @param vd.vars A character vector specifying the variables of interest to perform variance
#' decomposition.
#' @param assay A character string indicating the assay in the Seurat object to use for the
#' variance decomposition. Default is "GeneCooc".
#' @param cores An integer value specifying the number of CPU cores to use for parallel
#' computation. Default is 1.
#'
#' @return A data frame containing the variance decomposition results.
#'
#' @export
VarDecompose <- function(object, vd.vars, assay="GeneCooc", cores=1) {
  ## fetch data
  data <- t(GetAssayData(object, assay="GeneCooc", slot = "data"))
  meta.data <- object@meta.data
  genes <- colnames(data)
  if (!all(vd.vars %in% colnames(meta.data))) {
    vd.vars.404 <- setdiff(vd.vars, colnames(meta.data))
    stop(paste("vd.vars:", vd.vars.404, "is(are) not found in 'meta.data'"))
  }
  total.cores <- parallel::detectCores()
  cores <- min(cores, total.cores)
  ## prepare data
  vd.vars.str <- sapply(vd.vars, function(xx) sprintf("(1|%s)", xx))
  modelFormulaStr <- paste("expression ~", paste(vd.vars.str, collapse = "+"))
  data.use <- cbind(data[, genes], meta.data)
  ## run VD
  vd.res <- do.call(rbind, parallel::mclapply(genes, function(genename) {
    data.model <- data.use[, c(vd.vars, genename)]
    colnames(data.model) <- c(vd.vars, "expression")
    tryCatch({
      model <- suppressWarnings(lme4::lmer(stats::as.formula(modelFormulaStr), data = data.model, REML = TRUE, verbose = FALSE))
      results <- as.data.frame(lme4::VarCorr(model))
      rownames(results) <- results$grp
      results <- results[c(vd.vars, "Residual"), ]
      frac.var <- results$vcov / sum(results$vcov)

      res.tmp <- c("OK", frac.var)
    },
    error = function(e) {
      print(e)
      res.tmp <- c("FAIL", rep(-1, length(vd.vars)+1))
    })
    names(res.tmp) <- c("status", vd.vars, "residual")
    as.data.frame(as.list(res.tmp)) # return
  }, mc.cores = cores))
  rownames(vd.res) <- genes
  vd.res %<>% as.data.frame()
  vd.res <- vd.res %>% dplyr::mutate(module.name = rownames(vd.res), .before=1)
  for (i in 3:ncol(vd.res)) {
    vd.res[[i]] %<>% as.numeric()
  }
  return(vd.res)
}
