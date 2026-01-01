#' Generalized Pareto Distribution P-value Approximation
#'
#' @description Approximate p-value using GPD for extreme value statistics.
#' This is an internal function used by coptest and coptest.p.
#'
#' @param stat Numeric. The test statistic value.
#' @param null Numeric vector. The null distribution from permutation.
#' @param initial_p Numeric. Initial quantile threshold (default 0.8).
#'
#' @return Numeric. Approximated p-value.
#'
#' @keywords internal
#' @export
gpd_pval_approx <- function(stat, null, initial_p = 0.8) {
  if (length(null) == 0) {
    warning("Empty null distribution provided")
    return(NA_real_)
  }

  p <- initial_p
  q <- quantile(null, probs = p)
  d0 <- null[null >= q]
  d1 <- d0 - q

  # Handle edge case where d1 has too few observations

if (length(d1) < 3) {
    return(sum(null >= stat) / length(null))
  }

  gpd_test <- gPdtest::gpd.test(d1, J = 1000)

  # Adjust threshold if GPD fit is poor
  while (gpd_test$p.values[2, 1] <= 0.05 && p <= 0.99) {
    p <- p + 0.01
    q <- quantile(null, probs = p)
    d0 <- null[null >= q]
    d1 <- d0 - q

    if (length(d1) < 3) {
      return(sum(null >= stat) / length(null))
    }

    gpd_test <- gPdtest::gpd.test(d1)
  }

  gpd_fit <- gPdtest::gpd.fit(d1, method = "amle")
  k <- gpd_fit[1, 1]
  a <- gpd_fit[2, 1]

  # Handle edge cases in GPD calculation
  if (a <= 0 || is.na(k) || is.na(a)) {
    return(sum(null >= stat) / length(null))
  }

  f <- 1 - (1 + k * (stat - q) / a)^(-1 / k)
  pval <- (1 - f) * length(d0) / length(null)

  # Ensure p-value is in valid range
  pval <- max(0, min(1, pval))
  names(pval) <- NULL

  return(pval)
}


#' Validate Matrix Input for Copula Tests
#'
#' @description Validate that inputs are valid numeric matrices with compatible dimensions.
#'
#' @param x1 First input matrix.
#' @param x2 Second input matrix.
#' @param require_colnames Logical. Whether column names are required (default FALSE).
#'
#' @return NULL (invisible). Throws error if validation fails.
#'
#' @keywords internal
#' @export
validate_matrix_input <- function(x1, x2, require_colnames = FALSE) {
  # Check for NULL
  if (is.null(x1) || is.null(x2)) {
    stop("Inputs x1 and x2 must be numeric matrices, not NULL.", call. = FALSE)
  }

  # Check matrix type
  if (!is.matrix(x1) || !is.matrix(x2)) {
    stop("Inputs x1 and x2 must be numeric matrices.", call. = FALSE)
  }

  # Check numeric type
  if (!is.numeric(x1) || !is.numeric(x2)) {
    stop("Inputs x1 and x2 must contain numeric values.", call. = FALSE)
  }

  # Check column dimensions
  if (ncol(x1) != ncol(x2)) {
    stop(
      sprintf(
        "x1 and x2 must have the same number of columns. Got %d and %d.",
        ncol(x1), ncol(x2)
      ),
      call. = FALSE
    )
  }

  # Check minimum dimensions
  if (ncol(x1) < 2) {
    stop("Matrices must have at least 2 columns for copula testing.", call. = FALSE)
  }

  if (nrow(x1) < 3 || nrow(x2) < 3) {
    stop("Matrices must have at least 3 rows (samples) for copula testing.", call. = FALSE)
  }

  # Check for NA/NaN/Inf
  if (any(!is.finite(x1)) || any(!is.finite(x2))) {
    stop("Matrices contain NA, NaN, or Inf values. Please remove or impute these values.", call. = FALSE)
  }

  # Check column names if required
  if (require_colnames) {
    if (is.null(colnames(x1)) || is.null(colnames(x2))) {
      stop("Column names are required but missing.", call. = FALSE)
    }

    if (!all(colnames(x1) == colnames(x2))) {
      stop("Column names must match between x1 and x2.", call. = FALSE)
    }
  }

  invisible(NULL)
}
