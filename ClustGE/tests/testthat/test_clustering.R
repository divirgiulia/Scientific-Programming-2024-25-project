library(testthat)
library(SummarizedExperiment)
library(Biobase)
library(igraph)

# Create example data
set.seed(123)
mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(mat) <- paste0("Gene", 1:10)
colnames(mat) <- paste0("Sample", 1:10)
se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
exset <- Biobase::ExpressionSet(assayData = mat)

test_that("Hierarchical clustering returns hclust object", {
  hc <- hierarchical_clust(se)
  expect_s3_class(hc, "hclust")

  hc2 <- hierarchical_clust(exset, target = "genes")
  expect_s3_class(hc2, "hclust")
})

test_that("Hierarchical clustering with annotation adds cluster IDs", {
  result <- hierarchical_clust(se, annotate = TRUE, k = 3)
  expect_true(is.list(result))
  expect_s4_class(result$object, "SummarizedExperiment")
  expect_true("Cluster" %in% colnames(SummarizedExperiment::colData(result$object)))
})

test_that("K-means clustering returns kmeans object", {
  km <- kmeans_clust(se, centers = 3)
  expect_s3_class(km, "kmeans")

  km2 <- kmeans_clust(exset, centers = 2, target = "genes")
  expect_s3_class(km2, "kmeans")
})

test_that("K-means clustering with annotation adds cluster IDs", {
  result <- kmeans_clust(se, centers = 3, annotate = TRUE)
  expect_true(is.list(result))
  expect_s4_class(result$object, "SummarizedExperiment")
  expect_true("Cluster" %in% colnames(SummarizedExperiment::colData(result$object)))
})

test_that("Graph clustering returns igraph communities object", {
  comm <- graph_based_clust(se)
  expect_true(inherits(comm, "communities") || inherits(comm, "igraph_clusters"))

  comm2 <- graph_based_clust(exset, target = "genes")
  expect_true(inherits(comm2, "communities") || inherits(comm2, "igraph_clusters"))
})

test_that("Graph-based clustering with annotation adds cluster IDs", {
  result <- graph_based_clust(se, annotate = TRUE)
  expect_true(is.list(result))
  expect_s4_class(result$object, "SummarizedExperiment")
  expect_true("Cluster" %in% colnames(SummarizedExperiment::colData(result$object)))
})

test_that("clust_annotation works for different objects", {
  hc <- hclust(dist(t(mat)))
  se_annot <- clust_annotation(se, hc, k = 2)
  expect_true("Cluster" %in% colnames(SummarizedExperiment::colData(se_annot)))

  km <- kmeans(mat, centers = 2)
  exset_annot <- clust_annotation(exset, km, target = "samples")
  expect_true("Cluster" %in% colnames(Biobase::pData(exset_annot)))
})
