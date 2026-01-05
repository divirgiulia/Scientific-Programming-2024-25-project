test_that("get_expr_matrix extracts correctly", {
  mat <- matrix(1:6, nrow = 2)
  rownames(mat) <- paste0("gene", 1:nrow(mat))
  colnames(mat) <- paste0("sample", 1:ncol(mat))

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  eset <- Biobase::ExpressionSet(mat) # use shortcut, requires dimnames

  # SummarizedExperiment
  se_mat <- get_expr_matrix(se)
  expect_true(is.matrix(se_mat))
  expect_type(se_mat, "double")
  expect_equal(unname(se_mat), unname(mat))

  # ExpressionSet
  eset_mat <- get_expr_matrix(eset)
  expect_true(is.matrix(eset_mat))
  expect_type(eset_mat, "double")
  expect_equal(unname(eset_mat), unname(mat))

  # Wrong object
  expect_error(get_expr_matrix(list()))
})


test_that("set_expr_matrix replaces correctly", {
  mat <- matrix(1:6, nrow = 2)
  rownames(mat) <- paste0("gene", 1:nrow(mat))
  colnames(mat) <- paste0("sample", 1:ncol(mat))

  new_mat <- matrix(10:15, nrow = 2)
  rownames(new_mat) <- rownames(mat)
  colnames(new_mat) <- colnames(mat)

  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
  eset <- Biobase::ExpressionSet(mat)

  # SummarizedExperiment
  se2 <- set_expr_matrix(se, new_mat)
  expect_equal(unname(get_expr_matrix(se2)), unname(new_mat))

  # ExpressionSet
  eset2 <- set_expr_matrix(eset, new_mat)
  expect_equal(unname(get_expr_matrix(eset2)), unname(new_mat))

  # Wrong object
  expect_error(set_expr_matrix(list(), new_mat))
})
