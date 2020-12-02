#' Pairwise two sample test for empirical copula.
#' @description Autometed procedure to test the difference of copula with pairwise variables.
#' @param x1 Numeric matrix. Samples in row and variables in column
#' @param x2 Numeric matrix. The same with x1.
#' @param nperm The number of permutation.
#' @param approx Logical. If set `approx=TRUE`, p-value will be approximated using generalized Parato distribution. Otherwise, no approximation of p-value.
#' @param silent Logical. If FALSE, don't diplay progression message.
#' @details The difference of coptest is that coptest.p compare bivariate dependency structure rather than full joint dependencies as the coptest. Automatically generate all pairs of bivariate copula and compare them between the two conditions.
#' @return List of two components:
#' \itemize{
#'   \item tbl - data.frame consisting from varname(gene/protein names), stats(test statistic),p(p-value), p.adj(adjusted p-value).
#'   \item perm.out - list for information of randomization procedure. Each element is a variable pair. pval(p-value), stat(test statistics), p.approx(approximated p-value if approx=TRUE)
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
coptest.p = function(x1,x2,nperm=100,approx = TRUE,silent = TRUE){

  if(dim(x1)[2]!=dim(x2)[2]){
    stop("The number of columns is different between the conditions. Set the same dimensions.")
  }
  if(is.null(colnames(x1)) | is.null(colnames(x2))){
    stop("The column names are null")
  }
  if(!all(colnames(x1)==colnames(x2))){
    stop("The column names are diffrent between the conditions.")
  }

  comb = comb_n(ncol(x1),2)
  nc = ncol(comb)
  print(paste("There are",nc,"pairs variables."))

  #comb = comb_n(ncol(x),2)
  result_perm = vector("list",nc)
  names(result_perm) = apply(comb,2,function(ii)paste(colnames(x1)[ii[1]],colnames(x1)[ii[2]],sep="|"))

  for(i in 1:ncol(comb)){
    cidx = comb[,i]
    dice = coptest(x1[,cidx],x2[,cidx],nperm = nperm,approx = approx)
    pval = dice$pval
    stat = dice$stat
    null = dice$null
    result_perm[[i]]$pval = pval
    result_perm[[i]]$stat = stat
    result_perm[[i]]$null = null
    if(!silent){
      cat("Test for",i,"th /",nc,"pair of the variables was finished \n")
    }
  }
  cat("Two sample copula test was finished!\n")


  if(approx){
    cat("Approximating of p-value...\n")

    for(i in seq_along(result_perm)){
      stat = result_perm[[i]]$stat
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
        pval = (1 - f)*length(d0)/length(null)
        names(pval) = NULL
      }
      result_perm[[i]]$p.approx = pval
      cat("Approximating",i," th/",length(result_perm),"pairs of variable was finished!\n")
    }
    cat("Approximation was finished!!\n")

  }
  stats = sapply(result_perm,function(x)x$stat)
  pvals = sapply(result_perm,function(x)x$pval)
  qvals = p.adjust(pvals,method="BH")
  df = data.frame(varname = names(result_perm),stat = stats,p = pvals, p.adj = qvals,stringsAsFactors = FALSE)
  rownames(df) = NULL
  output = list(tbl = df, perm.out = result_perm)
}

