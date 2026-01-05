library(testthat)
library(SummarizedExperiment)
library(Biobase)
library(ggplot2)
library(pheatmap)

# Create test SummarizedExperiment
set.seed(123)
mat <- matrix(rnorm(50), nrow = 10, ncol = 5) # 10 genes x 5 samples
colnames(mat) <- paste0("Sample", 1:5)
rownames(mat) <- paste0("Gene", 1:10)
se <- SummarizedExperiment(assays = list(counts = mat))

# Add cluster annotation for samples and genes
colData(se)$Cluster <- factor(sample(1:2, ncol(se), replace = TRUE)) # 5 samples
rowData(se)$Cluster <- factor(sample(1:3, nrow(se), replace = TRUE)) # 10 genes

# Create test ExpressionSet
exset <- ExpressionSet(assayData = mat)
pData(exset)$Cluster <- factor(sample(1:2, ncol(exprs(exset)), replace = TRUE))
fData(exset)$Cluster <- factor(sample(1:3, nrow(exprs(exset)), replace = TRUE))


# heatmap
test_that("plot_clustered_heatmap works with SummarizedExperiment", {
  # Create mock expression matrix
  mat <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(mat) <- paste0("Gene", 1:3)
  colnames(mat) <- paste0("Sample", 1:4)

  # Create cluster annotations
  col_data <- DataFrame(Cluster = c("A", "A", "B", "B"))
  row_data <- DataFrame(Cluster = c("X", "Y", "X"))

  # Create SummarizedExperiment object
  se <- SummarizedExperiment(
    assays = list(counts = mat),
    colData = col_data,
    rowData = row_data
  )

  # Expect function runs without error
  expect_silent({
    p <- plot_clustered_heatmap(se, assay_name = "counts", target = "samples", cluster_col = "Cluster")
  })

  # The output is a pheatmap object
  expect_s3_class(p, "pheatmap")

  # Also test gene annotation
  expect_silent({
    p2 <- plot_clustered_heatmap(se, assay_name = "counts", target = "genes", cluster_col = "Cluster")
  })
  expect_s3_class(p2, "pheatmap")
})




# Create test SummarizedExperiment
set.seed(123)
mat <- matrix(rnorm(50), nrow = 10, ncol = 5) # 10 genes x 5 samples
colnames(mat) <- paste0("Sample", 1:5)
rownames(mat) <- paste0("Gene", 1:10)
se <- SummarizedExperiment(assays = list(counts = mat))
colData(se)$Cluster <- factor(sample(1:2, ncol(se), replace = TRUE))

# Create test ExpressionSet
exset <- ExpressionSet(assayData = mat)
pData(exset)$Cluster <- factor(sample(1:2, ncol(exprs(exset)), replace = TRUE))


# test pca
test_that("plot_pca returns a ggplot object", {
  # Without color_by
  expect_s3_class(plot_pca(se), "ggplot")
  expect_s3_class(plot_pca(exset), "ggplot")

  # With color_by
  expect_s3_class(plot_pca(se, color_by = "Cluster"), "ggplot")
  expect_s3_class(plot_pca(exset, color_by = "Cluster"), "ggplot")
})

test_that("plot_pca handles invalid color_by", {
  # Should not fail if color_by does not exist
  expect_s3_class(plot_pca(se, color_by = "Nonexistent"), "ggplot")
  expect_s3_class(plot_pca(exset, color_by = "Nonexistent"), "ggplot")
})
