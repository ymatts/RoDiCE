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
#'   \item stat - test statistics based on Cramér–von Mises test statistic
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
coptest = function(x1,x2,nperm=100,approx = TRUE){

  result = rcpp_coptest(mat1 = x1,mat2 = x2,nperm = nperm)
  pval = result$pval
  stat = result$stat
  null = result$null_stat

  if(approx){
    if(sum(null >= stat) <= 10){

      p = 0.8
      q = quantile(null,probs = p)
      d0 = null[null >= q]
      d1 = d0 - q
      gpd_test = gPdtest::gpd.test(d1,J = 1000)

      while(gpd_test$p.values[2,1] <= 0.05 & p <= 1){
        p = p + 0.01
        q = quantile(null,probs = p)
        d0 = null[null >= q]
        d1 = d0 - q
        gpd_test = gPdtest::gpd.test(d1)
      }

      gpd_fit = gPdtest::gpd.fit(d1,method = "amle")
      k = gpd_fit[1,1]
      a = gpd_fit[2,1]
      f = 1 - (1 + k * (stat - q) / a)^(-1 / k)
      #f = 1 - (1 - k * d1 / a)^(1 / k)
      pval = (1 - f)*length(d0)/length(null)
      names(pval) = NULL
    }
  }

  list(pval=pval,stat=stat,null=null)
}

