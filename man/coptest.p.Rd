% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coptest_p.R
\name{coptest.p}
\alias{coptest.p}
\title{Pairwise two sample test for empirical copula.}
\usage{
coptest.p(x1, x2, nperm = 100, approx = TRUE, silent = TRUE)
}
\arguments{
\item{x1}{Numeric matrix. Samples in row and variables in column}

\item{x2}{Numeric matrix. The same with x1.}

\item{nperm}{The number of permutation.}

\item{approx}{Logical. If set `approx=TRUE`, p-value will be approximated using generalized Parato distribution. Otherwise, no approximation of p-value.}

\item{silent}{Logical. If FALSE, don't diplay progression message.}
}
\value{
List of two components:
\itemize{
  \item tbl - data.frame consisting from varname(gene/protein names), stats(test statistic),p(p-value), p.adj(adjusted p-value).
  \item perm.out - list for information of randomization procedure. Each element is a variable pair. pval(p-value), stat(test statistics), p.approx(approximated p-value if approx=TRUE)
  \item stat - test statistics based on Cramér–von Mises test statistic
  \item null - null distribution from randomization procedure
}
}
\description{
Autometed procedure to test the difference of copula with pairwise variables.
}
\details{
The difference of coptest is that coptest.p compare bivariate dependency structure rather than full joint dependencies as the coptest. Automatically generate all pairs of bivariate copula and compare them between the two conditions.
}
\examples{
data(ccrcc.pbaf) # example data from clear renal cell carcinoma(clerk et al.2019)
data(corum.hsp.pbaf)
tumor = ccrcc.pbaf$tumor # 110 samples and 10 proteins from PBAF complex
normal = ccrcc.pbaf$normal # 84 samples and 10 proteins from PBAF complex

# perform copula test for pairwise variables.
result = coptest.p(tumor,normal,nperm=100,approx=TRUE)
result$tbl

}
\references{
{
Yusuke MATSUI et al.(2020) RoDiCE: Robust differential protein co-expression analysis for cancer complexome (submitted).

Clerk DJ et al.(2019) Integrated Proteogenomic Characterization of Clear Cell Renal Cell Carcinoma.Cell;179(4),964-983 e931.
}
}
\author{
Yusuke MATSUI
}
