#' Read and Convert Expression Data to SummarizedExperiment
#'
#' Reads gene expression data from a CSV file, ExpressionSet, or SummarizedExperiment
#' and returns a SummarizedExperiment object for downstream analysis.
#'
#' @param file Path to CSV file (rows = genes, columns = samples).
#' @param exset Optional: an ExpressionSet or SummarizedExperiment object.
#' @param expname Name of the assay slot for SummarizedExperiment input (default = "counts").
#' @param sample.annotation Optional data.frame with sample metadata (rows = samples).
#'
#' @return A SummarizedExperiment object.
#' @importFrom Biobase exprs
#' @importFrom SummarizedExperiment assay SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#' @importFrom utils read.csv
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' mat <- matrix(rpois(20, lambda = 10), nrow = 5)
#' rownames(mat) <- paste0("Gene", 1:5)
#' colnames(mat) <- paste0("Sample", 1:4)
#' write.csv(mat, tmp, row.names = TRUE)
#' se1 <- read_expression(file = tmp)
#' se1
#'
#' pdata <- data.frame(SampleGroup = c("A", "A", "B", "B"), row.names = colnames(mat))
#' eset <- Biobase::ExpressionSet(
#'   assayData = mat,
#'   phenoData = new("AnnotatedDataFrame", data = pdata)
#' )
#' se2 <- read_expression(exset = eset)
#' se2
#'
#' se0 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
#' se3 <- read_expression(exset = se0, expname = "counts")
#' se3
#'
#' sample_ann <- data.frame(
#'   Group = c("Ctrl", "Ctrl", "Treat", "Treat"),
#'   row.names = colnames(mat)
#' )
#' se4 <- read_expression(file = tmp, sample.annotation = sample_ann)
#' SummarizedExperiment::colData(se4)
#' @export
read_expression <- function(file = NULL, exset = NULL, expname = "counts", sample.annotation = NULL) {
  if (!is.null(exset)) {
    if (inherits(exset, "ExpressionSet")) {
      mat <- Biobase::exprs(exset)
    } else if (inherits(exset, "SummarizedExperiment")) {
      mat <- SummarizedExperiment::assay(exset, expname)
    } else {
      stop("exset must be an ExpressionSet or SummarizedExperiment")
    }
  } else if (!is.null(file)) {
    mat <- as.matrix(utils::read.csv(file, row.names = 1, check.names = FALSE))
  } else {
    stop("Provide either 'file' or 'exset'")
  }


  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = as.matrix(mat)))


  if (!is.null(sample.annotation)) {
    if (!all(colnames(mat) %in% rownames(sample.annotation))) {
      stop("Row names of sample.annotation must match column names of expression matrix")
    }
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(sample.annotation)
  }

  return(se)
}
