
#' Visualize differential protein complex co-expression network
#'
#' @description Create a plot of differential co-expressed proteins within a complex.
#' @param obj data.frame. Direct result object of `coptest.p`.
#' @param tbl data.frame. Alternative with `obj` when you manually prepare the table. This table must include variables three variables;`varname`(pairs of variable which is separated by `|`), `stat`(test statistic), `p`(p-value), and `p.adj`(adjusted p-value).
#' @param ref data.frame. Group information of the variables, e.g., protein complex. The table must include `group_id`(group identifier, e,g., protein complex identifier),`desc`(description of the group,e.g., protein complex name),and `varname`(variable name, e.g., proteins)
#' @param p Single numeric. Significance level for two sample test. The default is 0.05.
#' @param title character. Main title of the plot. if `ref` is provided, `title` is set to `desc`.
#'
#'@author Yusuke Matsui
#'
#'@references{
#'Yusuke MATSUI et al.(2020) RoDiCE: Robust differential protein co-expression analysis for cancer complexome (submitted).
#'}
#' @return Interactive network.
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
#' # without reference
#' netvis(result,title="PBAF complex")
#'
#' # with reference
#' data(corum.hsp.pbaf) # protein complex (PBAF complex; id=149) membership from CORUM database.
#' netvis(result,ref=corum.hsp.pbaf)
#'
#'@author Yusuke MATSUI
#'@references{
#'Yusuke MATSUI et al.(2020) RoDiCE: Robust differential protein co-expression analysis for cancer complexome (submitted).
#'}
#'

netvis = function(obj = NULL,tbl = NULL,ref = NULL, title = NULL,p = 0.05){
  if(is.null(obj) & is.null(tbl)){
    stop("set one of the arguments for `obj` or `tbl`.")
  }
  if(!is.null(obj) & is.null(tbl)){
    tbl = obj$tbl
  }
  if(is.null(obj) & !is.null(tbl)){
    tbl = tbl
  }

  if(is.null(ref)){

    varname = unique(unlist(strsplit(tbl$varname,split = "\\|")))
    id = seq_along(varname)
    var2id = data.frame(varname = varname,id = id)

    var_tbl = do.call(rbind,strsplit(tbl$varname,split = "\\|"))
    id_tbl = t(apply(var_tbl,1,function(x)id[match(x,varname)]))

    fullEl = cbind(id_tbl,var_tbl)
    colnames(fullEl) = c("from","to","from_label","to_label")

    graphDf = data.frame(fullEl)

    graphDf$group = 1 # 0: not expressed 1: expressed

    graphDf$pval = tbl$p[match(tbl$varname,paste(graphDf$from_label,graphDf$to_label,sep="|"))]
    graphDf$qval = tbl$p.adj[match(tbl$varname,paste(graphDf$from_label,graphDf$to_label,sep="|"))]
    graphDf$weight = -log10(graphDf$pval)
    graphDf$isSignifP = ifelse(graphDf$pval <= p,1,0)
    graphDf$isSignifQ = ifelse(graphDf$qval <= p,1,0)

    fullMember = unique(as.vector(fullEl[,3:4]))
    memberDf = data.frame(label = fullMember,group = 1)
    memberDf$id = id[match(memberDf$label,varname)]

    nodes = memberDf
    links = graphDf

    vis.nodes = nodes
    vis.links = links

    vis.nodes$size  = 10
    vis.nodes$shape  = "dot"
    vis.nodes$shadow = TRUE # Nodes will drop shadow
    vis.nodes$title  = nodes$label # Text on click
    vis.nodes$label  = nodes$label # Node labelF
    vis.nodes$borderWidth = 2 # Node border width
    vis.nodes$isExprs = ifelse(vis.nodes$group,"is.Expressed","isNOT.Expressed")

    colV = c("lightgrey","gold")
    vis.nodes$color.background = colV[vis.nodes$group+1]
    vis.nodes$color.border = "black"
    vis.nodes$color.highlight.background = "orange"
    vis.nodes$color.highlight.border = "darkred"

    colE = c('rgba(0,0,0,0)',"blue","tomato")
    # vis.links = vis.links[vis.links$group==1,]
    vis.links$width = 3*vis.links$weight*vis.links$isSignifQ # line width
    vis.links$color = colE[vis.links$group + vis.links$isSignifQ + 1]    # line color
    vis.links$arrows = "middle" # arrows: 'from', 'to', or 'middle'
    vis.links$smooth = FALSE    # should the edges be curved?
    vis.links$shadow = FALSE    # edge shadow
    vis.links$arrows = NA
    vis.links$length = 500

    # https://datastorm-open.github.io/visNetwork/legend.html
    lnode = data.frame(color = colV,label=c("Not Expressed","Expressed"))
    ledge = data.frame(color = colE[-1],label = c("Not significant","Significant"))

    visNetwork::visNetwork(vis.nodes, vis.links) %>%
      visInteraction(hover = TRUE,hideEdgesOnDrag = TRUE,selectConnectedEdges = TRUE) %>%
      visLegend(addEdges = ledge,addNodes = lnode,useGroups = FALSE) %>%
      visExport(type="pdf",label = "save") %>% addExport()

  }else{

    if(!all(c("group_id","desc","varname")==colnames(ref))){
      stop("Check the column name of the argument `ref`. You must set `group_id`,`desc`,`varname`. Follow the example code.")
    }

    comb = comb_n(nrow(ref),2)
    ref$var_id = as.numeric(factor(ref$varname))
    desc = ref$desc[1]
    group_id = ref$group_id[1]

    fullEl = matrix(NA,nrow=ncol(comb),ncol=4)
    colnames(fullEl) = c("from","to","from_label","to_label")

    for(i in 1:ncol(comb)){
      idx1 = comb[1,i]
      idx2 = comb[2,i]
      x1 = ref$var_id[idx1]
      x2 = ref$var_id[idx2]
      y1 = ref$varname[idx1]
      y2 = ref$varname[idx2]
      fullEl[i,] = c(x1,x2,y1,y2)
    }

    graphDf = data.frame(fullEl)
    graphDf$group = 0 # 0: not expressed 1: expressed

    varname = unique(unlist(strsplit(tbl$varname,split = "\\|")))
    id = seq_along(varname)
    var2id = data.frame(varname = varname,id = id)

    var_tbl = do.call(rbind,strsplit(tbl$varname,split = "\\|"))
    id_tbl = t(apply(var_tbl,1,function(x)id[match(x,varname)]))

    exprsIdx = rep(NA,nrow(var_tbl))
    for(i in 1:nrow(var_tbl)){
      r = var_tbl[i,]
      exprsIdx[i] = which(apply(fullEl[,3:4],1,function(x)all(x%in%r)))
    }

    graphDf$group[exprsIdx] = 1
    graphDf$pval = Inf
    graphDf$qval = Inf
    graphDf$pval[exprsIdx] = tbl$p
    graphDf$qval[exprsIdx] = tbl$p.adj

    graphDf$weight = -log10(graphDf$pval)
    graphDf$isSignifP = ifelse(graphDf$pval <= p,1,0)
    graphDf$isSignifQ = ifelse(graphDf$qval <= p,1,0)

    fullMember = unique(as.vector(fullEl[,3:4]))
    exprsMember = unique(as.vector(var_tbl))
    isExprsMember = fullMember %in% exprsMember
    memberDf = data.frame(label = fullMember,group = isExprsMember)
    memberDf$id = ref$var_id[match(memberDf$label,ref$varname)]

    nodes = memberDf
    links = graphDf

    vis.nodes = nodes
    vis.links = links

    vis.nodes$size  = 10
    vis.nodes$shape  = "dot"
    vis.nodes$shadow = TRUE # Nodes will drop shadow
    vis.nodes$title  = nodes$label # Text on click
    vis.nodes$label  = nodes$label # Node label
    vis.nodes$borderWidth = 2 # Node border width
    vis.nodes$isExprs = ifelse(vis.nodes$group,"is.Expressed","isNOT.Expressed")

    colV = c("lightgrey","gold")
    vis.nodes$color.background = colV[vis.nodes$group+1]
    vis.nodes$color.border = "black"
    vis.nodes$color.highlight.background = "orange"
    vis.nodes$color.highlight.border = "darkred"

    colE = c('rgba(0,0,0,0)',"blue","tomato")
    vis.links$width = 3*vis.links$weight*vis.links$isSignifQ # line width
    vis.links$color = colE[vis.links$group + vis.links$isSignifQ + 1]    # line color
    vis.links$arrows = "middle" # arrows: 'from', 'to', or 'middle'
    vis.links$smooth = FALSE    # should the edges be curved?
    vis.links$shadow = FALSE    # edge shadow
    vis.links$arrows = NA
    vis.links$length = 500

    lnode = data.frame(color = colV,label=c("Not Expressed","Expressed"))
    ledge = data.frame(color = colE[-1],label = c("Not significant","Significant"))

    visNetwork::visNetwork(vis.nodes, vis.links, main = desc) %>%
      visInteraction(hover = TRUE,hideEdgesOnDrag = TRUE,selectConnectedEdges = TRUE) %>%
      visLegend(addEdges = ledge,addNodes = lnode,useGroups = FALSE) %>%
      visExport(type="pdf",name = paste(group_id,desc,sep="_"),label = "save") %>% addExport()
  }

}
