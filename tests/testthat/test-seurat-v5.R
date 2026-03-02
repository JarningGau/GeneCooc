capture_warnings <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

make_test_seurat_v5 <- function(n_genes = 300, n_cells = 80) {
  set.seed(123)
  counts <- Matrix::Matrix(
    matrix(rpois(n_genes * n_cells, lambda = 2), nrow = n_genes, ncol = n_cells),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("cell", seq_len(n_cells))
  Seurat::CreateSeuratObject(counts = counts)
}

test_that("CalGeneRankings works on Seurat v5 object without slot deprecation warnings", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  suppressWarnings(skip_if_not_installed("CelliD"))

  obj <- make_test_seurat_v5()
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  res <- capture_warnings(
    CalGeneRankings(
      obj,
      min.expr.cells = 5,
      nfeatures = 200,
      ndim.mca = 10,
      use.variable.features = TRUE
    )
  )
  obj2 <- res$value
  expect_false(any(grepl("GeneCooc package", res$warnings, fixed = TRUE)))
  expect_true("GeneCooc" %in% names(obj2@misc))
  expect_true("features.use" %in% names(obj2@misc$GeneCooc))
  expect_true("gene.rankings" %in% names(obj2@misc$GeneCooc))
})

test_that("VarDecompose supports Seurat v5 assay access", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  skip_if_not_installed("lme4")

  obj <- make_test_seurat_v5(n_genes = 220, n_cells = 90)
  obj$batch <- sample(c("b1", "b2"), size = ncol(obj), replace = TRUE)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  layer_data <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
  suppressWarnings({
    obj[["GeneCooc"]] <- Seurat::CreateAssayObject(data = layer_data[1:200, , drop = FALSE])
  })

  res <- capture_warnings(VarDecompose(obj, vd.vars = "batch"))
  vd <- res$value
  expect_false(any(grepl("GeneCooc package", res$warnings, fixed = TRUE)))
  expect_true(is.data.frame(vd))
  expect_true("module.name" %in% colnames(vd))
  expect_true("batch" %in% colnames(vd))
})
