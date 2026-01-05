#' Plot Distribution of Gene Expression Values
#'
#' Creates a histogram of expression values from the counts assay.
#'
#' @param se A SummarizedExperiment or ExpressionSet object.
#' @param bins Number of bins for the histogram (default = 50).
#'
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom SummarizedExperiment assay
#' @importFrom Biobase exprs
#' @importFrom utils read.csv
#' @export
#' @examples
#' library(SummarizedExperiment)
#'
#' set.seed(1)
#' mat <- matrix(rpois(50, lambda = 10), nrow = 10)
#' rownames(mat) <- paste0("gene", 1:10)
#' colnames(mat) <- paste0("sample", 1:5)
#'
#' se <- SummarizedExperiment(assays = list(counts = mat))
#'
#' plot_expression_distribution(se, bins = 20)
plot_expression_distribution <- function(se, bins = 50) {
  if (inherits(se, "SummarizedExperiment")) {
    exprs_mat <- SummarizedExperiment::assay(se, "counts")
  } else if (inherits(se, "ExpressionSet")) {
    exprs_mat <- Biobase::exprs(se)
  } else {
    stop("Input must be either a SummarizedExperiment or ExpressionSet object")
  }

  df <- data.frame(value = as.vector(exprs_mat))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = "steelblue", color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Distribution of Expression Values",
      x = "Expression",
      y = "Frequency"
    )
}



#' Boxplot of Sample Expression Distributions
#'
#' Creates a boxplot showing the distribution of expression values for each sample.
#'
#' @param se A SummarizedExperiment or ExpressionSet object.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @importFrom SummarizedExperiment assay
#' @importFrom Biobase exprs
#' @importFrom rlang .data
#' @export
#' @examples
#' library(SummarizedExperiment)
#' set.seed(1)
#' mat <- matrix(rpois(50, lambda = 10), nrow = 10)
#' rownames(mat) <- paste0("gene", 1:10)
#' colnames(mat) <- paste0("sample", 1:5)
#' se <- SummarizedExperiment(assays = list(counts = mat))
#' plot_sample_boxplot(se)
plot_sample_boxplot <- function(se) {
  if (inherits(se, "SummarizedExperiment")) {
    exprs_mat <- SummarizedExperiment::assay(se, "counts")
  } else if (inherits(se, "ExpressionSet")) {
    exprs_mat <- Biobase::exprs(se)
  } else {
    stop("Input must be either a SummarizedExperiment or ExpressionSet object")
  }

  df <- as.data.frame(exprs_mat)
  df_long <- reshape2::melt(df, variable.name = "Sample", value.name = "Expression")

  ggplot2::ggplot(df_long, ggplot2::aes(x = .data$Sample, y = .data$Expression)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(
      title = "Expression Distribution per Sample",
      x = "Sample",
      y = "Expression"
    )
}
