
#' Visualize dependence patterns
#' @description This function visulaize dependence structure within a group using pairwise scatterplot of original scale and copula-transformed scale.
#' @param obj data.frame. Direct result object of `coptest.p`.
#' @param tbl data.frame. Alternative with `obj` when you manually prepare the table. This table must include variables three variables;`varname`(pairs of variable which is separated by `|`), `stat`(test statistic), `p`(p-value), and `p.adj`(adjusted p-value).
#' @param exprs data.frame. Matrix with variables in the rows and samples in columns.
#' @param grp character vector. The vector should be the same length with the number of rows of the `exprs` object.
#' @param p Single numeric. Significance level for two sample test. The default is 0.05.
#' @param title character. Main title of the plot. if `ref` is provided, `title` is set to `desc`.
#'
#' @details This plot can be used for investigating differential dependence structure in the copula-transformed space compared with the original scale. Pairwise co-expression patterns are visualzed as scatter plot based on ggpairs from GGally pacakge. There are two types of scatter plots in the upper and lower diagonal, respectively; lower one represents scatter plots with original scale, which indicates standard correaltion structure; second one represents the scatter plots copula-transformed scale, i.e., rank-based correlation structure. Smoothing is conducted with cubic spline.
#' @return See ggpairs\{GGally\} function.
#' @export
#'
#' @examples
#' data(ccrcc.pbaf) # example data from clear renal cell carcinoma(clerk et al.2019)
#' data(corum.hsp.pbaf)
#' tumor = ccrcc.pbaf$tumor # 110 samples and 10 proteins from PBAF complex
#' normal = ccrcc.pbaf$normal # 84 samples and 10 proteins from PBAF complex
#'
#' #perform copula test for pairwise variables.
#' result = coptest.p(tumor,normal,nperm=100,approx=TRUE)
#' result$tbl
#'
#' exprs = rbind(tumor,normal)
#' grp = c(rep(1,nrow(tumor)),rep(2,nrow(normal)))
#' coexvis(obj = result,exprs = exprs,grp = grp,p = 0.05, title = "PBAF complex")
#'
#'
#'@author Yusuke MATSUI
#'@references{
#'Yusuke MATSUI et al.(2020) RoDiCE: Robust differential protein co-expression analysis for cancer complexome (submitted).
#'}
#'

coexvis = function(obj = NULL,tbl = NULL,exprs=NULL,grp = NULL,p = 0.05,title = NULL){

  require(GGally)
  if(is.null(obj) & is.null(tbl)){
    stop("set one of the arguments for `obj` or `tbl`.")
  }
  if(!is.null(obj) & is.null(tbl)){
    tbl = obj$tbl
  }
  if(is.null(obj) & !is.null(tbl)){
    tbl = tbl
  }

  if(is.null(grp)){
    stop("set grp argument.\n")
  }

  if(length(grp)!=nrow(exprs)){
    stop("the length of grp is different from nrow(exprs).\n")
  }

  grp = factor(grp)
  tbl2 = tbl[tbl$p.adj <= p,]

  varname = unique(unlist(strsplit(tbl$varname,split = "\\|")))
  var_tbl = do.call(rbind,strsplit(tbl$varname,split = "\\|"))
  cidx = match(varname,colnames(exprs))

  exprs = exprs[,cidx,drop=F]
  colnames(exprs) = colnames(exprs)[cidx]

  df_exprs = data.frame(exprs,check.names = FALSE)
  df_exprs$grp = grp

  upperfun = function(data,mapping){

    grpid = unique(grp)
    for(i in seq_along(grpid)){
      tmp = data[data$grp==grpid[i],]
      mat = as.matrix(tmp[,-ncol(tmp)])
      pmat = rcpp_pobs(mat)
      tmp[,-ncol(tmp)] = pmat
      data[data$grp==grpid[i],] = tmp
    }

    ggplot(data = data, mapping = mapping) +
      geom_point(size=.5,alpha=0.5) +
      geom_smooth(method="gam",se=T)
  }

  pairs = ggpairs(df_exprs,aes(color=grp),
                  lower = list(continuous = wrap("smooth",size = .5, alpha = 0.3)),
                  upper = list(continuous = wrap(upperfun)),
                  diag = list(continuous = wrap("densityDiag", alpha = 0.5))
  ) + ggtitle(title)

  primary_var = colnames(df_exprs)[1:(ncol(df_exprs)-1)]
  pvar_pos = match(primary_var, pairs$xAxisLabels)
  plots = vector("list",pairs$ncol*length(pvar_pos))
  count = 1
  for(i in seq_along(pvar_pos)){
    for(j in 1:pairs$ncol){
      plots[[count]] = getPlot(pairs,i = pvar_pos[i],j = j)
      count = count + 1
    }
  }

  ggmatrix(
    plots,
    nrow = length(pvar_pos),
    ncol = pairs$ncol,
    xAxisLabels = pairs$xAxisLabels,
    yAxisLabels = primary_var
  )

}

