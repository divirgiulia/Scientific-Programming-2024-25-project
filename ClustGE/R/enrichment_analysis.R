#' Get genes by cluster
#'
#' @param se SummarizedExperiment or ExpressionSet object
#' @param target "genes" (clustering must be done on genes)
#' @return List of gene vectors, one per cluster
#' @export
#'
#' @examples
#' library(SummarizedExperiment)
#' mat <- matrix(rnorm(40), nrow = 10)
#' rownames(mat) <- paste0("gene", 1:10)
#' colnames(mat) <- paste0("sample", 1:4)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' SummarizedExperiment::rowData(se)$Cluster <- rep(1:2, each = 5)
#' get_genes_by_cluster(se, target = "genes")
get_genes_by_cluster <- function(se, target = c("genes", "samples")) {
  target <- match.arg(target)
  if (target != "genes") stop("GSEA is performed on genes, not samples.")

  cluster_col <- if (inherits(se, "SummarizedExperiment")) {
    SummarizedExperiment::rowData(se)$Cluster
  } else {
    Biobase::fData(se)$Cluster
  }

  gene_ids <- rownames(se)
  clusters <- unique(cluster_col)
  genes_by_cluster <- lapply(clusters, function(cl) gene_ids[cluster_col == cl])
  names(genes_by_cluster) <- paste0("Cluster_", clusters)
  genes_by_cluster
}


#' Perform GO enrichment analysis per gene cluster
#'
#' @param genes_by_cluster List of gene vectors (output of get_genes_by_cluster)
#' @param orgdb OrgDb object (e.g., org.Hs.eg.db)
#' @param ont Ontology for GO: "BP", "MF", "CC" (default "BP")
#' @param keyType Character; key type used in OrgDb (default = "SYMBOL")
#' @return List of enrichGO results per cluster
#' @export
#' @examples
#' library(org.Hs.eg.db)
#' clusters <- list(
#'   Cluster_1 = c("TP53", "EGFR", "BRCA1"),
#'   Cluster_2 = c("MTOR", "AKT1", "MAPK1")
#' )
#' results <- perform_cluster_go_enrichment(clusters, org.Hs.eg.db, ont = "BP")
perform_cluster_go_enrichment <- function(genes_by_cluster, orgdb, ont = "BP", keyType = "SYMBOL") {
  lapply(genes_by_cluster, function(genes) {
    if (length(genes) > 0) {
      clusterProfiler::enrichGO(
        gene = genes,
        OrgDb = orgdb,
        keyType = keyType,
        ont = ont,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05,
        readable = TRUE
      )
    } else {
      NULL
    }
  })
}


#' Plot top GO terms per cluster
#'
#' @param results Named list of enrichGO objects (output of perform_cluster_go_enrichment)
#' @param n Integer; number of top GO terms to show per cluster (default = 10)
#' @return Invisibly returns the input results list
#' @export
#' @examples
#' library(org.Hs.eg.db)
#' clusters <- list(
#'   Cluster_1 = c("TP53", "EGFR", "BRCA1"),
#'   Cluster_2 = c("MTOR", "AKT1", "MAPK1")
#' )
#' results <- perform_cluster_go_enrichment(clusters, org.Hs.eg.db)
#' plot_cluster_go_enrichment(results, n = 5)
plot_cluster_go_enrichment <- function(results, n = 10) {
  if (is.null(results) || length(results) == 0) {
    warning("No enrichment results provided.")
    return(invisible(NULL))
  }

  lapply(names(results), function(cl) {
    res <- results[[cl]]
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      print(enrichplot::dotplot(res, showCategory = n) + ggplot2::ggtitle(cl))
    }
  })

  invisible(results)
}


#' Perform KEGG enrichment analysis on gene clusters
#'
#' @param se SummarizedExperiment or ExpressionSet object with cluster annotations.
#' @param target "genes" — KEGG enrichment is performed on gene clusters.
#' @param organism KEGG organism code (e.g., "hsa" for human, "mmu" for mouse).
#' @param keyType Type of gene identifiers ("kegg", "ncbi-geneid", or "ENTREZID").
#' @param pAdjustMethod Method for multiple testing correction (default "BH").
#' @param pvalueCutoff P-value cutoff (default 0.05).
#' @param qvalueCutoff Q-value cutoff (default 0.2).
#' @return List of enrichKEGG results per cluster.
#' @export
#' @examples
#' library(SummarizedExperiment)
#' mat <- matrix(rpois(40, 10), nrow = 10)
#' rownames(mat) <- paste0("gene", 1:10)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' SummarizedExperiment::rowData(se)$Cluster <- rep(1:2, each = 5)
#' results <- perform_enrichKEGG(se, organism = "hsa", keyType = "ncbi-geneid")
perform_enrichKEGG <- function(se, target = "genes",
                               organism, keyType = "kegg",
                               pAdjustMethod = "BH", pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2) {
  if (target != "genes") stop("KEGG enrichment is performed on genes, not samples.")

  if (inherits(se, "SummarizedExperiment")) {
    clusters <- SummarizedExperiment::rowData(se)$Cluster
    gene_ids <- rownames(se)
  } else if (inherits(se, "ExpressionSet")) {
    clusters <- Biobase::fData(se)$Cluster
    gene_ids <- rownames(se)
  } else {
    stop("Object must be a SummarizedExperiment or ExpressionSet")
  }

  cluster_levels <- unique(clusters)
  enrichment_results <- list()

  for (cl in cluster_levels) {
    genes_in_cluster <- gene_ids[clusters == cl]
    if (length(genes_in_cluster) > 0) {
      enrichment_results[[paste0("Cluster_", cl)]] <- clusterProfiler::enrichKEGG(
        gene = genes_in_cluster,
        organism = organism,
        keyType = keyType,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff
      )
    }
  }

  enrichment_results
}
