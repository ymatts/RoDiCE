#' Pairwise two sample test for empirical copula.
#' @description Automated procedure to test the difference of copula with pairwise variables.
#' @param x1 Numeric matrix. Samples in row and variables in column
#' @param x2 Numeric matrix. The same with x1.
#' @param nperm The number of permutation.
#' @param approx Logical. If set `approx=TRUE`, p-value will be approximated using generalized Parato distribution. Otherwise, no approximation of p-value.
#' @param silent Logical. If TRUE, don't display progression message.
#' @details The difference of coptest is that coptest.p compare bivariate dependency structure rather than full joint dependencies as the coptest. Automatically generate all pairs of bivariate copula and compare them between the two conditions.
#' @return List of two components:
#' \itemize{
#'   \item tbl - data.frame consisting from varname(gene/protein names), stats(test statistic),p(p-value), p.adj(adjusted p-value).
#'   \item perm.out - list for information of randomization procedure. Each element is a variable pair. pval(p-value), stat(test statistics), p.approx(approximated p-value if approx=TRUE)
#' }
#'
#' @export
#'
#' @examples
#' data(ccrcc.pbaf) # example data from clear renal cell carcinoma(clerk et al.2019)
#' data(corum.hsp.pbaf)
#' tumor = ccrcc.pbaf$tumor # 110 samples and 10 proteins from PBAF complex
#' normal = ccrcc.pbaf$normal # 84 samples and 10 proteins from PBAF complex
#'
#' # perform copula test for pairwise variables.
#' result = coptest.p(tumor,normal,nperm=100,approx=TRUE)
#' result$tbl
#'
#'@author Yusuke MATSUI
#'@references{
#'Yusuke MATSUI et al.(2020) RoDiCE: Robust differential protein co-expression analysis for cancer complexome (submitted).
#'
#'Clerk DJ et al.(2019) Integrated Proteogenomic Characterization of Clear Cell Renal Cell Carcinoma.Cell;179(4),964-983 e931.
#'}
coptest.p <- function(x1, x2, nperm = 100, approx = TRUE, silent = TRUE) {
  # Input validation
  validate_matrix_input(x1, x2, require_colnames = TRUE)

  # Validate nperm
  if (!is.numeric(nperm) || nperm < 1) {
    stop("nperm must be a positive integer.", call. = FALSE)
  }
  nperm <- as.integer(nperm)

  # Generate all pairwise combinations
  comb <- comb_n(ncol(x1), 2)
  nc <- ncol(comb)

  if (!silent) {
    message(sprintf("Testing %d variable pairs.", nc))
  }

  # Initialize result list
  result_perm <- vector("list", nc)
  names(result_perm) <- apply(comb, 2, function(ii) {
    paste(colnames(x1)[ii[1]], colnames(x1)[ii[2]], sep = "|")
  })

  # Run pairwise copula tests
  for (i in seq_len(ncol(comb))) {
    cidx <- comb[, i]
    # Note: approx=FALSE here since we handle approximation separately below
    dice <- coptest(x1[, cidx], x2[, cidx], nperm = nperm, approx = FALSE)

    result_perm[[i]]$pval <- dice$pval
    result_perm[[i]]$stat <- dice$stat
    result_perm[[i]]$null <- dice$null

    if (!silent) {
      message(sprintf("Test %d/%d completed.", i, nc))
    }
  }

  if (!silent) {
    message("Two sample copula test completed.")
  }

  # Apply GPD approximation if requested
  if (approx) {
    if (!silent) {
      message("Approximating p-values...")
    }

    for (i in seq_along(result_perm)) {
      stat <- result_perm[[i]]$stat
      null <- result_perm[[i]]$null  # Bug fix: use correct null distribution

      if (sum(null >= stat) <= 10) {
        pval_approx <- gpd_pval_approx(stat, null)
        result_perm[[i]]$pval <- pval_approx
      }

      result_perm[[i]]$p.approx <- result_perm[[i]]$pval

      if (!silent) {
        message(sprintf("Approximation %d/%d completed.", i, length(result_perm)))
      }
    }

    if (!silent) {
      message("P-value approximation completed.")
    }
  }

  # Compile results into data frame
  stats <- sapply(result_perm, function(x) x$stat)
  pvals <- sapply(result_perm, function(x) x$pval)
  qvals <- p.adjust(pvals, method = "BH")

  df <- data.frame(
    varname = names(result_perm),
    stat = stats,
    p = pvals,
    p.adj = qvals,
    stringsAsFactors = FALSE
  )
  rownames(df) <- NULL

  list(tbl = df, perm.out = result_perm)
}
