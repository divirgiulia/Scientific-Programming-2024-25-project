#' Remove/flag outlier genes based on mean expression
#'
#' Flags genes whose mean expression lies outside 1.5×IQR from the first or third quartile.
#' Optionally removes them.
#'
#' @param se A SummarizedExperiment or ExpressionSet object.
#' @param assay_name Character; assay to use for SummarizedExperiment (default "counts").
#' @param k Numeric; multiplier for IQR fence (default 1.5).
#' @param remove Logical; if TRUE, drop outlier genes (default FALSE).
#'
#' @return An object of the same type as input, with `Outlier` column in rowData/fData.
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom Biobase exprs fData
#' @export
#' @examples
#' library(SummarizedExperiment)
#' set.seed(123)
#' mat <- matrix(rpois(40, lambda = 10), nrow = 8)
#' rownames(mat) <- paste0("gene", 1:8)
#' colnames(mat) <- paste0("sample", 1:5)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' se_flagged <- remove_outlier_genes(se, assay_name = "counts", remove = FALSE)
#' head(SummarizedExperiment::rowData(se_flagged))
#' se_clean <- remove_outlier_genes(se, assay_name = "counts", remove = TRUE)
#' dim(se_clean)
remove_outlier_genes <- function(se, assay_name = "counts", k = 1.5, remove = FALSE) {
  if (inherits(se, "SummarizedExperiment")) {
    if (!assay_name %in% names(SummarizedExperiment::assays(se))) {
      stop(sprintf("Assay '%s' not found in SummarizedExperiment", assay_name))
    }
    exprs_mat <- SummarizedExperiment::assay(se, assay_name)
  } else if (inherits(se, "ExpressionSet")) {
    exprs_mat <- Biobase::exprs(se)
  } else {
    stop("Input must be SummarizedExperiment or ExpressionSet")
  }


  gene_means <- rowMeans(exprs_mat, na.rm = TRUE)


  q1 <- stats::quantile(gene_means, 0.25, na.rm = TRUE)
  q3 <- stats::quantile(gene_means, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - k * iqr
  upper <- q3 + k * iqr

  outlier_flag <- (gene_means < lower) | (gene_means > upper)


  if (inherits(se, "SummarizedExperiment")) {
    SummarizedExperiment::rowData(se)$Outlier <- outlier_flag
  } else {
    Biobase::fData(se)$Outlier <- outlier_flag
  }


  if (remove) {
    keep <- !outlier_flag
    if (inherits(se, "SummarizedExperiment")) {
      se <- se[keep, ]
    } else {
      se <- se[keep, ]
    }
  }

  se
}
