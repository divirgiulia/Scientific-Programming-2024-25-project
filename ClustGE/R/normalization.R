#' Normalize Expression Data
#'
#' Applies normalization to gene expression data.
#'
#' @param se A SummarizedExperiment or ExpressionSet object.
#' @param method Normalization method: "log" (log2(x+1)) or "quantile".
#' @param standardize Logical; if TRUE, scale genes to mean 0, sd 1.
#'
#' @return A normalized object of the same type as input.
#' @importFrom SummarizedExperiment assay
#' @importFrom stats prcomp sd
#' @importFrom Biobase exprs
#' @importFrom preprocessCore normalize.quantiles
#' @export
#' @examples
#' library(SummarizedExperiment)
#' mat <- matrix(rpois(40, lambda = 10), nrow = 10)
#' rownames(mat) <- paste0("gene", 1:10)
#' colnames(mat) <- paste0("sample", 1:4)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#'
#' se_norm_log <- normalize_expr(se, method = "log", standardize = TRUE)
#'
#' se_norm_q <- normalize_expr(se, method = "quantile", standardize = FALSE)
normalize_expr <- function(se, method = c("log", "quantile"), standardize = TRUE) {
  method <- match.arg(method)

  if (inherits(se, "SummarizedExperiment")) {
    exprs_mat <- SummarizedExperiment::assay(se)
  } else if (inherits(se, "ExpressionSet")) {
    exprs_mat <- Biobase::exprs(se)
  } else {
    stop("Input must be either a SummarizedExperiment or ExpressionSet object")
  }

  # normalization
  if (method == "log") {
    exprs_mat <- log2(exprs_mat + 1)
  } else if (method == "quantile") {
    exprs_mat <- preprocessCore::normalize.quantiles(exprs_mat)
    dimnames(exprs_mat) <- dimnames(if (inherits(se, "SummarizedExperiment")) {
      SummarizedExperiment::assay(se)
    } else {
      Biobase::exprs(se)
    })
  }


  if (standardize) {
    exprs_mat <- t(apply(exprs_mat, 1, function(x) {
      s <- sd(x)
      if (s == 0) rep(0, length(x)) else (x - mean(x)) / s
    }))
  }


  if (inherits(se, "SummarizedExperiment")) {
    SummarizedExperiment::assay(se) <- exprs_mat
  } else {
    Biobase::exprs(se) <- exprs_mat
  }
  return(se)
}
