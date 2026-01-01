context("coptest.p - Pairwise copula test")

test_that("coptest.p returns expected structure", {
  skip_on_cran()

  set.seed(123)
  x1 <- matrix(rnorm(50 * 4), nrow = 50, ncol = 4)
  x2 <- matrix(rnorm(40 * 4), nrow = 40, ncol = 4)
  colnames(x1) <- colnames(x2) <- c("A", "B", "C", "D")

  result <- coptest.p(x1, x2, nperm = 5, approx = FALSE, silent = TRUE)

  # Check return structure
  expect_true(is.list(result))
  expect_true(all(c("tbl", "perm.out") %in% names(result)))

  # Check tbl structure
  expect_true(is.data.frame(result$tbl))
  expect_true(all(c("varname", "stat", "p", "p.adj") %in% colnames(result$tbl)))

  # Number of pairs should be C(4,2) = 6
  expect_equal(nrow(result$tbl), 6)
})

test_that("coptest.p validates input dimensions", {
  set.seed(123)
  x1 <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
  x2 <- matrix(rnorm(40 * 4), nrow = 40, ncol = 4)
  colnames(x1) <- c("A", "B", "C")
  colnames(x2) <- c("A", "B", "C", "D")

  expect_error(coptest.p(x1, x2), "number of columns")
})

test_that("coptest.p validates column names", {
  set.seed(123)
  x1 <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
  x2 <- matrix(rnorm(40 * 3), nrow = 40, ncol = 3)

  # No column names
  expect_error(coptest.p(x1, x2), "column names are null")

  # Different column names
  colnames(x1) <- c("A", "B", "C")
  colnames(x2) <- c("X", "Y", "Z")
  expect_error(coptest.p(x1, x2), "column names are diffrent")
})
