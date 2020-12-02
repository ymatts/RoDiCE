#' Construct subunit-level expression data
#' @description Split whole protein expression into the subunits.
#' @param exprs data.frame. Protein expression data where variables(genes / proteins) in the columns and samples in the rows.
#' @param ref data.frame. Group information of the variables, e.g., protein complex. The table must include `group_id`(group identifier, e,g., protein complex identifier),`desc`(description of the group,e.g., protein complex name),and `varname`(variable name, e.g., proteins)
#'
#' @return list of which element corresponds to each complex.
#'
#' \itemize{
#'   \item comp.exprs - data.frame. Variables (genes/proteins) in the columns and samples in the rows.
#'   \item comp.names - character vector. Complex names
#'   \item comp.ids - data.frame. Raw data from CORUM database "allComplexes.txt" from http://mips.helmholtz-muenchen.de/corum/#download.
#' }
#'
#' @export
#' @author Yusuke Matsui
#'
#'@references{
#'Yusuke MATSUI et al.(2020) RoDiCE: Robust differential protein co-expression analysis for cancer complexome (submitted).
#'}

#' @examples
#' data(ccrcc.subset)
#' data(corum.subset)
#' tumor = ccrcc.subset$tumor
#' normal = ccrcc.subset$normal
#' tumorByComplex = group.exprs(tumor,corum.subset)
#' normalByComplex = group.exprs(normal,corum.subset)
#'
group.exprs = function(exprs,ref = NULL){

  unique_id = unique(ref$group_id)
  unique_name = ref$desc[match(unique_id,ref$group_id)]
  exprs.list = vector("list",length(unique_id))
  #names(exprs.list) = paste(unique_id,unique_name,sep = "|")

  for(i in seq_along(exprs.list)){
    id = unique_id[i]
    select_ref = ref[ref$group_id==id,]
    gene_names = select_ref$varname
    exprs.list[[i]] = exprs[,colnames(exprs) %in% gene_names,drop = FALSE]
  }

  output = list(comp.exprs = exprs.list, comp.names = unique_name, comp.ids = unique_id)
  output
}
