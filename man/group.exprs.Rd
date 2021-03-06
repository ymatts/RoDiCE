% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_exprs.R
\name{group.exprs}
\alias{group.exprs}
\title{Construct subunit-level expression data}
\usage{
group.exprs(exprs, ref = NULL)
}
\arguments{
\item{exprs}{data.frame. Protein expression data where variables(genes / proteins) in the columns and samples in the rows.}

\item{ref}{data.frame. Group information of the variables, e.g., protein complex. The table must include `group_id`(group identifier, e,g., protein complex identifier),`desc`(description of the group,e.g., protein complex name),and `varname`(variable name, e.g., proteins)}
}
\value{
list of which element corresponds to each complex.

\itemize{
  \item comp.exprs - data.frame. Variables (genes/proteins) in the columns and samples in the rows.
  \item comp.names - character vector. Complex names
  \item comp.ids - data.frame. Raw data from CORUM database "allComplexes.txt" from http://mips.helmholtz-muenchen.de/corum/#download.
}
}
\description{
Split whole protein expression into the subunits.
}
\examples{
data(ccrcc.subset)
data(corum.subset)
tumor = ccrcc.subset$tumor
normal = ccrcc.subset$normal
tumorByComplex = group.exprs(tumor,corum.subset)
normalByComplex = group.exprs(normal,corum.subset)

}
\references{
{
Yusuke MATSUI et al.(2020) RoDiCE: Robust differential protein co-expression analysis for cancer complexome (submitted).
}
}
\author{
Yusuke Matsui
}
