#' Two sample test for empirical copula.
#'
#' @param x1 Numeric matrix. Samples in row and variables in column
#' @param x2 Numeric matrix. The same with x1.
#' @param nperm The number of permutation.
#' @param approx Logical. If `approx=TRUE`, p-value will be approximated using generalized Parato distribution. Otherwise, no approximation of p-value.
#'
#' @return List of three components:
#' \itemize{
#'   \item pval - p-value
#'   \item stat - test statistics based on Cramer-von Mises test statistic
#'   \item null - null distribution from randomization procedure
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
#' # multivariate copula test(more than three variables)
#' result = coptest(tumor,normal,nperm=100,approx=FALSE)
#' result$pval
#'
#'@author Yusuke MATSUI
#'@references{
#'Yusuke MATSUI et al.(2020) RoDiCE: Robust differential protein co-expression analysis for cancer complexome (submitted).
#'
#'Clerk DJ et al.(2019) Integrated Proteogenomic Characterization of Clear Cell Renal Cell Carcinoma.Cell;179(4),964-983 e931.
#'}
coptest <- function(x1, x2, nperm = 100, approx = TRUE) {
  # Validate input
  validate_matrix_input(x1, x2, require_colnames = FALSE)

  # Validate nperm
  if (!is.numeric(nperm) || nperm < 1) {
    stop("nperm must be a positive integer.", call. = FALSE)
  }
  nperm <- as.integer(nperm)

  # Run C++ copula test
  result <- rcpp_coptest(mat1 = x1, mat2 = x2, nperm = nperm)
  pval <- result$pval
  stat <- result$stat
  null <- result$null_stat

  # Apply GPD approximation for extreme p-values if requested
  if (approx && sum(null >= stat) <= 10) {
    pval <- gpd_pval_approx(stat, null)
  }

  list(pval = pval, stat = stat, null = null)
}
