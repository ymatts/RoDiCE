context("coptest - Two sample copula test")

test_that("coptest returns expected structure", {
  skip_on_cran()

  # Create test data
  set.seed(123)
  x1 <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)
  x2 <- matrix(rnorm(80 * 3), nrow = 80, ncol = 3)
  colnames(x1) <- colnames(x2) <- c("A", "B", "C")

  result <- coptest(x1, x2, nperm = 10, approx = FALSE)

  # Check return structure

expect_true(is.list(result))
  expect_true(all(c("pval", "stat", "null") %in% names(result)))

  # Check value types
  expect_true(is.numeric(result$pval))
  expect_true(is.numeric(result$stat))
  expect_true(is.numeric(result$null))

  # Check value ranges
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_true(result$stat >= 0)
  expect_equal(length(result$null), 10)
})

test_that("coptest handles different sample sizes", {
  skip_on_cran()

  set.seed(456)
  x1 <- matrix(rnorm(50 * 2), nrow = 50, ncol = 2)
  x2 <- matrix(rnorm(30 * 2), nrow = 30, ncol = 2)

  result <- coptest(x1, x2, nperm = 5, approx = FALSE)

  expect_true(is.list(result))
  expect_true(result$pval >= 0 && result$pval <= 1)
})

test_that("coptest with approx=TRUE works", {
  skip_on_cran()

  set.seed(789)
  x1 <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)
  x2 <- matrix(rnorm(80 * 3), nrow = 80, ncol = 3)

  result <- coptest(x1, x2, nperm = 50, approx = TRUE)

  expect_true(is.list(result))
  expect_true(result$pval >= 0 && result$pval <= 1)
})
