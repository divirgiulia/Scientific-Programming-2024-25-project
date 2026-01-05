## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(ClustGE)
library(SummarizedExperiment)

## ----reading expression-------------------------------------------------------
library(ClustGE)
se <- read_expression(file = system.file("extdata", "gene_expr_example.csv", package = "ClustGE"))

## ----plots--------------------------------------------------------------------
plot_expression_distribution(se)
plot_sample_boxplot(se)

## ----normalization------------------------------------------------------------
library(SummarizedExperiment)
se_norm <- normalize_expr(se, method = "log", standardize = FALSE)
assay(se_norm, "counts")
plot_sample_boxplot(se_norm)

## ----outliers-----------------------------------------------------------------
se_clean <- remove_outlier_genes(se_norm, assay_name = "counts", remove = TRUE)

## ----kmeans annotated---------------------------------------------------------
km <- kmeans_clust(se_norm, centers = 5, target = "genes", annotate = TRUE)

## ----kmeans-------------------------------------------------------------------
km <- kmeans_clust(se_norm, centers = 5, target = "genes", annotate = FALSE)

## ----hclustering--------------------------------------------------------------
hc_annot <- hierarchical_clust(se_clean, target = "genes", annotate = TRUE, k = 4)
plot_clustered_heatmap(hc_annot$object, target = "genes", cluster_col = "Cluster")

## ----pca----------------------------------------------------------------------
hc_annot <- hierarchical_clust(se_clean, target = "samples", annotate = TRUE, k = 3)
plot_pca(hc_annot$object, color_by = "Cluster")

## ----dendrograms--------------------------------------------------------------
hc_annot <- hierarchical_clust(se_clean, target = "samples", annotate = TRUE, k = 3)
mat <- SummarizedExperiment::assay(hc_annot$object, "counts")
hc <- hclust(dist(t(mat)))
plot_dendrogram(hc, k = 3)

## ----clusters&genes-----------------------------------------------------------
hc_annot <- hierarchical_clust(se_clean, target = "genes", annotate = TRUE, k = 4)
genes_by_cluster <- get_genes_by_cluster(hc_annot$object, target = "genes")

## ----GO, message=FALSE, warning=FALSE-----------------------------------------
library(org.Hs.eg.db)
go_results <- perform_cluster_go_enrichment(
  genes_by_cluster,
  orgdb = org.Hs.eg.db,
  ont = "BP"
)
plot_cluster_go_enrichment(go_results, n = 10)

## ----KEGG, message=FALSE, warning=FALSE---------------------------------------
kegg_results <- perform_enrichKEGG(
  se_clean,
  target = "genes",
  organism = "hsa", # human
  keyType = "ncbi-geneid"
)

## ----session-info-------------------------------------------------------------
sessionInfo()

