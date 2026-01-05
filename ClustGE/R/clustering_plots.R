#' Plot heatmap with cluster annotation
#'
#' @param se SummarizedExperiment or ExpressionSet object
#' @param assay_name Name of assay
#' @param target "samples" or "genes" for annotation
#' @param cluster_col Name of column storing cluster IDs
#' @param scale Scale rows ("row") or not ("none")
#' @return pheatmap object
#' @importFrom pheatmap pheatmap
#' @importFrom Biobase featureNames sampleNames
#' @importFrom Biobase pData
#' @export
#' @examples
#' library(SummarizedExperiment)
#' mat <- matrix(rnorm(100), 10, 10)
#' colData <- DataFrame(Cluster = rep(1:2, each = 5))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = mat),
#'   colData = colData
#' )
#' plot_clustered_heatmap(se, target = "samples", cluster_col = "Cluster")
plot_clustered_heatmap <- function(se, assay_name = "counts", target = c("samples", "genes"),
                                   cluster_col = "Cluster", scale = "row") {
  target <- match.arg(target)
  mat <- get_expr_matrix(se, assay_name)

  if (is.null(rownames(mat))) {
    rownames(mat) <- if (inherits(se, "ExpressionSet")) Biobase::featureNames(se) else paste0("Gene", seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- if (inherits(se, "ExpressionSet")) Biobase::sampleNames(se) else paste0("Sample", seq_len(ncol(mat)))
  }

  annotation <- NULL
  if (inherits(se, "SummarizedExperiment")) {
    if (target == "samples") {
      annotation_vec <- colData(se)[[cluster_col]]
      if (length(annotation_vec) != ncol(mat)) stop("Cluster vector length does not match number of samples")
      annotation <- data.frame(Cluster = annotation_vec)
      rownames(annotation) <- colnames(mat)
    } else {
      annotation_vec <- rowData(se)[[cluster_col]]
      if (length(annotation_vec) != nrow(mat)) stop("Cluster vector length does not match number of genes")
      annotation <- data.frame(Cluster = annotation_vec)
      rownames(annotation) <- rownames(mat)
    }
  } else if (inherits(se, "ExpressionSet")) {
    if (target == "samples") {
      annotation_vec <- pData(se)[[cluster_col]]
      if (length(annotation_vec) != ncol(mat)) stop("Cluster vector length does not match number of samples")
      annotation <- data.frame(Cluster = annotation_vec)
      rownames(annotation) <- colnames(mat)
    } else {
      annotation_vec <- fData(se)[[cluster_col]]
      if (length(annotation_vec) != nrow(mat)) stop("Cluster vector length does not match number of genes")
      annotation <- data.frame(Cluster = annotation_vec)
      rownames(annotation) <- rownames(mat)
    }
  }

  pheatmap::pheatmap(mat,
    scale = scale,
    cluster_rows = TRUE, cluster_cols = TRUE,
    annotation_col = if (target == "samples") annotation else NULL,
    annotation_row = if (target == "genes") annotation else NULL
  )
}




#' PCA Plot of Samples #' #' Performs Principal Component Analysis on the counts assay and plots the first two principal components.
#'
#' @param se A SummarizedExperiment or ExpressionSet object.
#' @param color_by Column name in sample metadata (colData/pData) to color points (optional).
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom SummarizedExperiment assay colData
#' @importFrom Biobase exprs pData
#' @importFrom stats prcomp
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' library(SummarizedExperiment)
#' library(S4Vectors)
#' mat <- matrix(rnorm(100), 10, 10)
#' colData <- S4Vectors::DataFrame(Group = rep(c("A", "B"), each = 5))
#' se <- SummarizedExperiment(assays = list(counts = mat), colData = colData)
#' plot_pca(se, color_by = "Group")
plot_pca <- function(se, color_by = NULL) {
  if (inherits(se, "SummarizedExperiment")) {
    exprs_mat <- SummarizedExperiment::assay(se, "counts")
  } else if (inherits(se, "ExpressionSet")) {
    exprs_mat <- Biobase::exprs(se)
  } else {
    stop("Input must be either a SummarizedExperiment or ExpressionSet object")
  }

  pca <- stats::prcomp(t(exprs_mat), scale. = TRUE)
  pca_df <- as.data.frame(pca$x[, 1:2])
  colnames(pca_df) <- c("PC1", "PC2")

  if (!is.null(color_by)) {
    if (inherits(se, "SummarizedExperiment")) {
      pca_df[[color_by]] <- SummarizedExperiment::colData(se)[[color_by]]
    } else if (inherits(se, "ExpressionSet")) {
      pca_df[[color_by]] <- Biobase::pData(se)[[color_by]]
    }
  }

  ggplot2::ggplot(pca_df, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = if (!is.null(color_by)) .data[[color_by]] else NULL)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "PCA of Samples",
      x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)")
    )
}

#' Plot dendrogram of hierarchical clustering
#'
#' @param hc hclust object
#' @param k Optional, cut tree into k clusters
#' @return Dendrogram plot
#' @importFrom stats as.dendrogram
#' @importFrom ggdendro ggdendrogram
#' @import ggplot2
#' @export
#' @examples
#' mat <- matrix(rnorm(25), 5, 5)
#' hc <- hclust(dist(mat))
#' plot_dendrogram(hc, k = 2)
plot_dendrogram <- function(hc, k = NULL) {
  dend <- as.dendrogram(hc)
  if (!is.null(k)) dend <- dendextend::color_branches(dend, k = k)
  ggdendro::ggdendrogram(dend, rotate = FALSE) +
    theme_minimal() +
    labs(title = "Hierarchical Clustering Dendrogram")
}
