library(testthat)
library(SummarizedExperiment)
library(Biobase)



test_that("get_genes_by_cluster extracts genes correctly (SummarizedExperiment)", {
  mat <- matrix(rnorm(20), nrow = 5)
  rownames(mat) <- paste0("Gene", 1:5)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat),
    rowData = data.frame(Cluster = c(1, 1, 2, 2, 3))
  )

  genes_by_cluster <- get_genes_by_cluster(se, target = "genes")
  expect_type(genes_by_cluster, "list")
  expect_true(all(grepl("Cluster_", names(genes_by_cluster))))
  expect_equal(length(genes_by_cluster), 3)
})




test_that("get_genes_by_cluster extracts genes correctly (ExpressionSet)", {
  mat <- matrix(rnorm(20), nrow = 5)
  rownames(mat) <- paste0("Gene", 1:5)

  # Ensure featureData rownames match assayData rownames
  fdat <- data.frame(
    Cluster = c(1, 1, 2, 2, 3),
    row.names = rownames(mat) # must match exactly
  )

  eset <- ExpressionSet(
    assayData = mat,
    featureData = AnnotatedDataFrame(fdat)
  )

  genes_by_cluster <- get_genes_by_cluster(eset, target = "genes")

  expect_type(genes_by_cluster, "list")
  expect_true(all(grepl("Cluster_", names(genes_by_cluster))))
  expect_equal(length(genes_by_cluster), 3)

  # Check cluster contents
  expect_equal(genes_by_cluster$Cluster_1, c("Gene1", "Gene2"))
  expect_equal(genes_by_cluster$Cluster_2, c("Gene3", "Gene4"))
  expect_equal(genes_by_cluster$Cluster_3, "Gene5")
})




test_that("get_genes_by_cluster errors for target = 'samples'", {
  mat <- matrix(rnorm(20), nrow = 5)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat),
    rowData = data.frame(Cluster = 1:5)
  )
  expect_error(get_genes_by_cluster(se, target = "samples"))
})




# GO enrichment test (mocked small gene set + org.Hs.eg.db if available)
test_that("perform_cluster_go_enrichment returns a list of enrichGO objects or NULL", {
  skip_if_not_installed("org.Hs.eg.db")

  # load the OrgDb object explicitly
  orgdb <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))

  genes_by_cluster <- list(
    Cluster_1 = c("TP53", "EGFR"),
    Cluster_2 = c("BRCA1", "BRCA2")
  )

  results <- perform_cluster_go_enrichment(genes_by_cluster, orgdb)

  expect_type(results, "list")
  expect_true(all(names(results) %in% names(genes_by_cluster)))
  expect_true(all(vapply(results, function(x) is.null(x) || inherits(x, "enrichResult"), logical(1))))
})




# test Plot_cluster_go_enrichment
test_that("plot_cluster_go_enrichment runs and returns list invisibly", {
  skip_if_not_installed("org.Hs.eg.db")
  orgdb <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))

  genes_by_cluster <- list(Cluster_1 = c("TP53", "EGFR"))
  results <- perform_cluster_go_enrichment(genes_by_cluster, orgdb)

  expect_silent({
    plots <- plot_cluster_go_enrichment(results, n = 5)
    expect_type(plots, "list")
    expect_true(all(names(plots) %in% names(results)))
  })
})




# KEGG enrichment test (uses dummy gene IDs)
test_that("perform_enrichKEGG returns a list of enrichKEGG objects", {
  skip_if_not_installed("clusterProfiler")
  # KEGG expects ENTREZ IDs usually
  mat <- matrix(rnorm(20), nrow = 5)
  rownames(mat) <- c("7157", "1956", "2064", "673", "5290") # TP53, EGFR, ERBB2, BRAF, PIK3CA
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat),
    rowData = data.frame(Cluster = c(1, 1, 2, 2, 3))
  )

  results <- perform_enrichKEGG(se, organism = "hsa", keyType = "ncbi-geneid")
  expect_type(results, "list")
  expect_true(all(grepl("Cluster_", names(results))))
})
