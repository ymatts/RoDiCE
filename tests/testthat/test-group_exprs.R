context("group.exprs - Split expression data by complex")

test_that("group.exprs returns expected structure", {
  # Create test data
  exprs <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  colnames(exprs) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")

  ref <- data.frame(
    group_id = c("C1", "C1", "C2", "C2", "C2"),
    desc = c("Complex1", "Complex1", "Complex2", "Complex2", "Complex2"),
    varname = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5"),
    stringsAsFactors = FALSE
  )

  result <- group.exprs(exprs, ref)

  # Check return structure
  expect_true(is.list(result))
  expect_true(all(c("comp.exprs", "comp.names", "comp.ids") %in% names(result)))

  # Check number of complexes
  expect_equal(length(result$comp.exprs), 2)
  expect_equal(length(result$comp.names), 2)
  expect_equal(length(result$comp.ids), 2)

  # Check complex membership
  expect_equal(ncol(result$comp.exprs[[1]]), 2)  # Complex1 has 2 genes
  expect_equal(ncol(result$comp.exprs[[2]]), 3)  # Complex2 has 3 genes
})

test_that("group.exprs handles missing genes gracefully", {
  exprs <- matrix(rnorm(20 * 3), nrow = 20, ncol = 3)
  colnames(exprs) <- c("Gene1", "Gene2", "Gene3")

  ref <- data.frame(
    group_id = c("C1", "C1"),
    desc = c("Complex1", "Complex1"),
    varname = c("Gene1", "GeneX"),  # GeneX not in exprs
    stringsAsFactors = FALSE
  )

  result <- group.exprs(exprs, ref)

  # Should only include Gene1
  expect_equal(ncol(result$comp.exprs[[1]]), 1)
})
