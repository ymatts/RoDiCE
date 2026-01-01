context("Utility functions")

test_that("gpd_pval_approx returns valid p-value", {
  skip_on_cran()

  set.seed(123)
  null <- rnorm(1000, mean = 0, sd = 1)
  stat <- 2.5

  result <- gpd_pval_approx(stat, null)

  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("gpd_pval_approx handles edge cases", {
  skip_on_cran()

  set.seed(456)
  null <- rnorm(1000)
  stat <- max(null) + 1  # stat larger than all null values

  result <- gpd_pval_approx(stat, null)

  expect_true(is.numeric(result))
  expect_true(result >= 0)
})

test_that("validate_matrix_input detects invalid input", {
  # NULL input
  expect_error(validate_matrix_input(NULL, NULL), "must be numeric matrices")

  # Non-matrix input
  expect_error(validate_matrix_input(c(1, 2, 3), c(4, 5, 6)), "must be numeric matrices")

  # Different number of columns
  x1 <- matrix(1:6, nrow = 2, ncol = 3)
  x2 <- matrix(1:8, nrow = 2, ncol = 4)
  expect_error(validate_matrix_input(x1, x2), "same number of columns")
})

test_that("validate_matrix_input accepts valid input", {
  x1 <- matrix(rnorm(20), nrow = 10, ncol = 2)
  x2 <- matrix(rnorm(16), nrow = 8, ncol = 2)

  expect_silent(validate_matrix_input(x1, x2))
})
