#' Extract expression matrix from object
#'
#' A helper function to extract the expression matrix from a
#' SummarizedExperiment or ExpressionSet object.
#'
#' @param se A SummarizedExperiment or ExpressionSet object.
#' @param assay_name Character. The name of the assay to use.
#'
#' @return A numeric matrix of expression values.
#' @importFrom SummarizedExperiment assay
#' @importFrom Biobase exprs
#' @export
#' @examples
#' mat <- matrix(rpois(12, lambda = 10), nrow = 3)
#' rownames(mat) <- paste0("gene", 1:3)
#' colnames(mat) <- paste0("sample", 1:4)
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = mat)
#' )
#' get_expr_matrix(se, assay_name = "counts")
get_expr_matrix <- function(se, assay_name = "counts") {
  if (inherits(se, "SummarizedExperiment")) {
    mat <- SummarizedExperiment::assay(se, assay_name)
  } else if (inherits(se, "ExpressionSet")) {
    mat <- Biobase::exprs(se)
  } else {
    stop("Object must be a SummarizedExperiment or ExpressionSet")
  }


  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  mat
}


#' Replce expression matrix in object
#'
#' @param se A SummarizedExperiment or ExpressionSet object.
#' @param exprs_mat A numeric matrix to replace the expression data.
#' @param assay_name Name of assay to replace.
#'
#' @return Updated object of the same type as input.
#' @export
#' @examples
#' mat <- matrix(rpois(12, lambda = 10), nrow = 3)
#' rownames(mat) <- paste0("gene", 1:3)
#' colnames(mat) <- paste0("sample", 1:4)
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = mat)
#' )
#' new_mat <- mat + 1
#' se <- set_expr_matrix(se, new_mat, assay_name = "counts")
#' get_expr_matrix(se, "counts")
set_expr_matrix <- function(se, exprs_mat, assay_name = "counts") {
  if (inherits(se, "SummarizedExperiment")) {
    SummarizedExperiment::assay(se, assay_name) <- exprs_mat
  } else if (inherits(se, "ExpressionSet")) {
    Biobase::exprs(se) <- exprs_mat
  } else {
    stop("Object must be a SummarizedExperiment or ExpressionSet")
  }
  se
}
