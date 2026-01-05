library(testthat)
library(SummarizedExperiment)
library(ClustGE)

context("Preprocessing functions")

# Test read_expression()
test_that("read_expression loads a CSV correctly", {
  # create temporary CSV
  tmpfile <- tempfile(fileext = ".csv")
  mat <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(mat) <- paste0("Sample", 1:4)
  rownames(mat) <- paste0("Gene", 1:3)
  write.csv(mat, tmpfile, row.names = TRUE)

  se <- read_expression(tmpfile)

  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(dim(se), dim(mat))
  expect_equal(rownames(se), rownames(mat))
  expect_equal(colnames(se), colnames(mat))
})


# Test normalize_expression()
test_that("normalize_expression scales data correctly", {
  mat <- matrix(1:9, nrow = 3)
  se <- SummarizedExperiment(assays = list(counts = mat))

  se_norm <- normalize_expr(se)

  expect_s4_class(se_norm, "SummarizedExperiment")
  expect_equal(dim(assay(se_norm)), dim(mat))
  # Check mean ~0, sd ~1 if standardization
  scaled <- assay(se_norm)
  expect_true(abs(mean(scaled)) < 1e-8)
  expect_true(all(abs(apply(scaled, 1, sd) - 1) < 1e-8))
})


# Test remove_outlier_genes()
test_that("remove_outlier_genes flags and optionally removes by IQR", {
  # 4 genes x 5 samples; one very high-mean gene
  mat <- matrix(
    c(
      1, 2, 3, 4, 5,
      2, 2, 2, 2, 2,
      3, 3, 3, 3, 3,
      1000, 1000, 1000, 1000, 1000
    ),
    nrow = 4, byrow = TRUE
  )
  se <- SummarizedExperiment(assays = list(counts = mat))

  # flag only
  se_flag <- remove_outlier_genes(se, k = 1.5, remove = FALSE)
  expect_s4_class(se_flag, "SummarizedExperiment")
  expect_true("Outlier" %in% colnames(rowData(se_flag)))
  expect_equal(sum(rowData(se_flag)$Outlier), 1)

  # remove
  se_removed <- remove_outlier_genes(se, k = 1.5, remove = TRUE)
  expect_equal(nrow(se_removed), 3)
})
