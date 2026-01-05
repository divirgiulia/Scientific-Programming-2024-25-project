#' Perform hierarchical clustering on expression data
#'
#' @param se SummarizedExperiment or ExpressionSet object.
#' @param method Distance method (default "euclidean").
#' @param clust_method Clustering method (default "complete").
#' @param assay_name Name of assay (SummarizedExperiment only).
#' @param target "samples" or "genes" — whether to cluster samples or genes.
#' @param annotate Logical; if TRUE, add cluster IDs to object metadata.
#' @param k Number of clusters (only used if annotate = TRUE and method is hclust).
#' @return hclust object, or list(object, clustering) if annotate = TRUE.
#' @importFrom stats dist hclust
#' @export
#' @examples
#' library(SummarizedExperiment)
#' mat <- matrix(rnorm(50), nrow = 10) # 10 genes x 5 samples
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' hc <- hierarchical_clust(se, target = "samples")
#' res <- hierarchical_clust(se, target = "samples", annotate = TRUE, k = 2)
#' names(res)
hierarchical_clust <- function(se, method = "euclidean", clust_method = "complete",
                               assay_name = "counts", target = c("samples", "genes"),
                               annotate = FALSE, k = NULL) {
  target <- match.arg(target)
  exprs_mat <- get_expr_matrix(se, assay_name)

  dist_mat <- if (target == "samples") {
    dist(t(exprs_mat), method = method)
  } else {
    dist(exprs_mat, method = method)
  }

  hc <- hclust(dist_mat, method = clust_method)

  if (annotate) {
    se <- clust_annotation(se, hc, k = k, target = target, assay_name = assay_name)
    return(list(object = se, clustering = hc))
  }

  return(hc)
}


#' Perform k-means clustering on expression data
#'
#' @param se SummarizedExperiment or ExpressionSet object.
#' @param centers Number of clusters.
#' @param assay_name Name of assay (SummarizedExperiment only).
#' @param target "samples" or "genes" — whether to cluster samples or genes.
#' @param nstart Number of random starts (default 10).
#' @param annotate Logical; if TRUE, add cluster IDs to object metadata.
#' @return kmeans object, or list(object, clustering) if annotate = TRUE.
#' @importFrom stats kmeans
#' @export
#' @examples
#' library(SummarizedExperiment)
#' mat <- matrix(rnorm(50), nrow = 10)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' km <- kmeans_clust(se, centers = 2, target = "samples")
#' res <- kmeans_clust(se, centers = 2, target = "samples", annotate = TRUE)
#' names(res)
kmeans_clust <- function(se, centers, assay_name = "counts", target = c("samples", "genes"),
                         nstart = 10, annotate = FALSE) {
  target <- match.arg(target)
  exprs_mat <- get_expr_matrix(se, assay_name)

  km <- if (target == "samples") {
    kmeans(t(exprs_mat), centers = centers, nstart = nstart)
  } else {
    kmeans(exprs_mat, centers = centers, nstart = nstart)
  }

  if (annotate) {
    se <- clust_annotation(se, km, target = target, assay_name = assay_name)
    return(list(object = se, clustering = km))
  }

  return(km)
}


#' Perform graph-based clustering (Louvain community detection) on expression data
#'
#' @param se SummarizedExperiment or ExpressionSet object.
#' @param k Number of nearest neighbors for graph (default 5).
#' @param assay_name Name of assay (SummarizedExperiment only).
#' @param target "samples" or "genes".
#' @param annotate Logical; if TRUE, add cluster IDs to object metadata.
#' @return igraph communities object, or list(object, clustering) if annotate = TRUE.
#' @importFrom igraph graph_from_edgelist cluster_louvain membership
#' @importFrom FNN get.knn
#' @export
#' @examples
#' library(SummarizedExperiment)
#' mat <- matrix(rnorm(50), nrow = 10)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' comm <- graph_based_clust(se, target = "samples", k = 3)
#' res <- graph_based_clust(se, target = "samples", k = 3, annotate = TRUE)
#' names(res)
graph_based_clust <- function(se, k = 5, assay_name = "counts",
                              target = c("samples", "genes"), annotate = FALSE) {
  target <- match.arg(target)
  exprs_mat <- get_expr_matrix(se, assay_name)

  mat <- if (target == "samples") t(exprs_mat) else exprs_mat
  nn <- FNN::get.knn(mat, k = k)
  edges <- cbind(rep(1:nrow(mat), each = k), as.vector(nn$nn.index))
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)

  comm <- igraph::cluster_louvain(g)

  if (annotate) {
    se <- clust_annotation(se, comm, target = target)
    return(list(object = se, clustering = comm))
  }

  return(comm)
}


#' Annotate samples or genes with cluster assignments
#'
#' @description
#' Adds cluster IDs from various clustering results to a SummarizedExperiment or ExpressionSet object.
#'
#' @param se SummarizedExperiment or ExpressionSet object.
#' @param clustering_result Result from hclust, kmeans, dbscan, PAM, or igraph community detection.
#' @param k Number of clusters to cut (only for hclust; ignored otherwise).
#' @param target "samples" or "genes" — whether clustering was done on samples or genes.
#' @param assay_name Name of assay (SummarizedExperiment only).
#' @param cluster_col Column name in metadata to use for annotation
#' @return Updated object with cluster assignments stored in metadata.
#' @importFrom stats cutree
#' @export
#' @examples
#' library(SummarizedExperiment)
#' mat <- matrix(rnorm(50), nrow = 10)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' hc <- hclust(dist(t(mat)))
#' se_annot <- clust_annotation(se, hc, k = 2, target = "samples")
#' SummarizedExperiment::colData(se_annot)
clust_annotation <- function(se, clustering_result, k = NULL,
                             target = c("samples", "genes"),
                             assay_name = "counts", cluster_col = "Cluster") {
  target <- match.arg(target)

  cluster_ids <- NULL

  if (inherits(clustering_result, "hclust")) {
    if (is.null(k)) stop("For hclust objects, please specify k (number of clusters)")
    cluster_ids <- stats::cutree(clustering_result, k = k)
  } else if (inherits(clustering_result, "kmeans")) {
    cluster_ids <- clustering_result$cluster
  } else if (inherits(clustering_result, "dbscan")) {
    cluster_ids <- clustering_result$cluster
    cluster_ids[cluster_ids == 0] <- NA
  } else if (inherits(clustering_result, "pam")) {
    cluster_ids <- clustering_result$clustering
  } else if (inherits(clustering_result, "communities") || inherits(clustering_result, "igraph_clusters")) {
    cluster_ids <- igraph::membership(clustering_result)
  } else {
    stop("Unsupported clustering result. Must be from hclust, kmeans, dbscan, PAM, or igraph communities.")
  }

  cluster_ids <- factor(cluster_ids, levels = sort(unique(cluster_ids)))

  if (inherits(se, "SummarizedExperiment")) {
    if (target == "samples") {
      SummarizedExperiment::colData(se)[[cluster_col]] <- cluster_ids
    } else {
      SummarizedExperiment::rowData(se)[[cluster_col]] <- cluster_ids
    }
  } else if (inherits(se, "ExpressionSet")) {
    if (target == "samples") {
      Biobase::pData(se)[[cluster_col]] <- cluster_ids
    } else {
      Biobase::fData(se)[[cluster_col]] <- cluster_ids
    }
  } else {
    stop("Object must be SummarizedExperiment or ExpressionSet")
  }

  return(se)
}
