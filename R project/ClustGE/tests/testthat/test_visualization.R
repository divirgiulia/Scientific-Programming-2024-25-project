# testing plot_expression_distribution() (data exploration)
test_that("plot_expression_distribution produces a ggplot", {
  # test SummarizedExperiment input
  mat <- matrix(1:10, nrow = 5)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))

  p <- plot_expression_distribution(se, bins = 5)

  expect_s3_class(p, "ggplot")
  expect_true(any(sapply(p$layers, function(l) inherits(l$geom, "Geom"))))

  # test ExpressionSet input
  eset <- Biobase::ExpressionSet(assayData = mat)

  p2 <- plot_expression_distribution(eset, bins = 5)

  expect_s3_class(p2, "ggplot")
  expect_true(any(sapply(p$layers, function(l) inherits(l$geom, "Geom"))))
})


# sample boxplot test
test_that("plot_sample_boxplot produces a ggplot with boxplot geom", {
  # toy expression matrix
  mat <- matrix(1:12, nrow = 3)
  rownames(mat) <- paste0("gene", 1:3)
  colnames(mat) <- paste0("sample", 1:4)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  eset <- Biobase::ExpressionSet(mat)

  # SummarizedExperiment input
  p1 <- plot_sample_boxplot(se)
  expect_s3_class(p1, "ggplot")
  expect_true(any(sapply(p1$layers, function(l) inherits(l$geom, "GeomBoxplot"))))

  # ExpressionSet input
  p2 <- plot_sample_boxplot(eset)
  expect_s3_class(p2, "ggplot")
  expect_true(any(sapply(p2$layers, function(l) inherits(l$geom, "GeomBoxplot"))))

  # Wrong input type should error
  expect_error(plot_sample_boxplot(list()))
})
