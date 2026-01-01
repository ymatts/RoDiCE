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

netvis <- function(obj = NULL, tbl = NULL, ref = NULL, title = NULL, p = 0.05) {
  # Input validation
  tbl <- .validate_netvis_input(obj, tbl)

  # Validate p threshold
  if (!is.numeric(p) || p <= 0 || p > 1) {
    stop("p must be a numeric value between 0 and 1.", call. = FALSE)
  }

  if (is.null(ref)) {
    .netvis_without_ref(tbl, title, p)
  } else {
    .netvis_with_ref(tbl, ref, p)
  }
}


#' Validate netvis input arguments
#' @keywords internal
.validate_netvis_input <- function(obj, tbl) {
  if (is.null(obj) && is.null(tbl)) {
    stop("Please provide either 'obj' (coptest.p result) or 'tbl' (data.frame).", call. = FALSE)
  }

  if (!is.null(obj) && is.null(tbl)) {
    if (!is.list(obj) || !"tbl" %in% names(obj)) {
      stop("'obj' must be a result from coptest.p containing a 'tbl' element.", call. = FALSE)
    }
    tbl <- obj$tbl
  }

  # Validate required columns
  required_cols <- c("varname", "stat", "p", "p.adj")
  missing_cols <- setdiff(required_cols, colnames(tbl))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in tbl: %s", paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  tbl
}


#' Style nodes for visualization
#' @keywords internal
.style_nodes <- function(nodes) {
  colV <- c("lightgrey", "gold")

  nodes$size <- 10
  nodes$shape <- "dot"
  nodes$shadow <- TRUE
  nodes$title <- nodes$label
  nodes$borderWidth <- 2
  nodes$isExprs <- ifelse(nodes$group, "is.Expressed", "isNOT.Expressed")
  nodes$color.background <- colV[nodes$group + 1]
  nodes$color.border <- "black"
  nodes$color.highlight.background <- "orange"
  nodes$color.highlight.border <- "darkred"

  nodes
}


#' Style links for visualization
#' @keywords internal
.style_links <- function(links) {
  colE <- c("rgba(0,0,0,0)", "blue", "tomato")

  links$width <- 3 * links$weight * links$isSignifQ
  links$color <- colE[links$group + links$isSignifQ + 1]
  links$smooth <- FALSE
  links$shadow <- FALSE
  links$arrows <- NA
  links$length <- 500

  links
}


#' Create legend data
#' @keywords internal
.create_legend <- function() {
  colV <- c("lightgrey", "gold")
  colE <- c("rgba(0,0,0,0)", "blue", "tomato")

  lnode <- data.frame(color = colV, label = c("Not Expressed", "Expressed"))
  ledge <- data.frame(color = colE[-1], label = c("Not significant", "Significant"))

  list(nodes = lnode, edges = ledge)
}


#' Network visualization without reference
#' @keywords internal
.netvis_without_ref <- function(tbl, title, p) {
  # Parse variable names
  varname <- unique(unlist(strsplit(tbl$varname, split = "\\|")))
  id <- seq_along(varname)

  var_tbl <- do.call(rbind, strsplit(tbl$varname, split = "\\|"))
  id_tbl <- t(apply(var_tbl, 1, function(x) id[match(x, varname)]))

  fullEl <- cbind(id_tbl, var_tbl)
  colnames(fullEl) <- c("from", "to", "from_label", "to_label")

  # Build graph data frame
  graphDf <- data.frame(fullEl, stringsAsFactors = FALSE)
  graphDf$group <- 1

  graphDf$pval <- tbl$p[match(tbl$varname, paste(graphDf$from_label, graphDf$to_label, sep = "|"))]
  graphDf$qval <- tbl$p.adj[match(tbl$varname, paste(graphDf$from_label, graphDf$to_label, sep = "|"))]
  graphDf$weight <- -log10(graphDf$pval)
  graphDf$isSignifP <- ifelse(graphDf$pval <= p, 1, 0)
  graphDf$isSignifQ <- ifelse(graphDf$qval <= p, 1, 0)

  # Build member data frame
  fullMember <- unique(as.vector(fullEl[, 3:4]))
  memberDf <- data.frame(label = fullMember, group = 1, stringsAsFactors = FALSE)
  memberDf$id <- id[match(memberDf$label, varname)]

  # Apply styling
  vis.nodes <- .style_nodes(memberDf)
  vis.links <- .style_links(graphDf)
  legend <- .create_legend()

  # Create visualization
  visNetwork::visNetwork(vis.nodes, vis.links, main = title) %>%
    visInteraction(hover = TRUE, hideEdgesOnDrag = TRUE, selectConnectedEdges = TRUE) %>%
    visLegend(addEdges = legend$edges, addNodes = legend$nodes, useGroups = FALSE) %>%
    visExport(type = "pdf", label = "save") %>%
    addExport()
}


#' Network visualization with reference
#' @keywords internal
.netvis_with_ref <- function(tbl, ref, p) {
  # Validate reference columns
  required_ref_cols <- c("group_id", "desc", "varname")
  if (!all(required_ref_cols %in% colnames(ref))) {
    stop(sprintf(
      "Reference must contain columns: %s. Found: %s",
      paste(required_ref_cols, collapse = ", "),
      paste(colnames(ref), collapse = ", ")
    ), call. = FALSE)
  }

  # Generate all possible pairs from reference
  comb <- comb_n(nrow(ref), 2)
  ref$var_id <- as.numeric(factor(ref$varname))
  desc <- ref$desc[1]
  group_id <- ref$group_id[1]

  fullEl <- matrix(NA, nrow = ncol(comb), ncol = 4)
  colnames(fullEl) <- c("from", "to", "from_label", "to_label")

  for (i in seq_len(ncol(comb))) {
    idx1 <- comb[1, i]
    idx2 <- comb[2, i]
    fullEl[i, ] <- c(ref$var_id[idx1], ref$var_id[idx2],
                     ref$varname[idx1], ref$varname[idx2])
  }

  # Build graph data frame
  graphDf <- data.frame(fullEl, stringsAsFactors = FALSE)
  graphDf$group <- 0

  # Match with test results
  varname <- unique(unlist(strsplit(tbl$varname, split = "\\|")))
  var_tbl <- do.call(rbind, strsplit(tbl$varname, split = "\\|"))

  exprsIdx <- sapply(seq_len(nrow(var_tbl)), function(i) {
    r <- var_tbl[i, ]
    which(apply(fullEl[, 3:4], 1, function(x) all(x %in% r)))
  })

  graphDf$group[exprsIdx] <- 1
  graphDf$pval <- Inf
  graphDf$qval <- Inf
  graphDf$pval[exprsIdx] <- tbl$p
  graphDf$qval[exprsIdx] <- tbl$p.adj

  graphDf$weight <- -log10(graphDf$pval)
  graphDf$isSignifP <- ifelse(graphDf$pval <= p, 1, 0)
  graphDf$isSignifQ <- ifelse(graphDf$qval <= p, 1, 0)

  # Build member data frame
  fullMember <- unique(as.vector(fullEl[, 3:4]))
  exprsMember <- unique(as.vector(var_tbl))
  isExprsMember <- fullMember %in% exprsMember
  memberDf <- data.frame(label = fullMember, group = isExprsMember, stringsAsFactors = FALSE)
  memberDf$id <- ref$var_id[match(memberDf$label, ref$varname)]

  # Apply styling
  vis.nodes <- .style_nodes(memberDf)
  vis.links <- .style_links(graphDf)
  legend <- .create_legend()

  # Create visualization
  visNetwork::visNetwork(vis.nodes, vis.links, main = desc) %>%
    visInteraction(hover = TRUE, hideEdgesOnDrag = TRUE, selectConnectedEdges = TRUE) %>%
    visLegend(addEdges = legend$edges, addNodes = legend$nodes, useGroups = FALSE) %>%
    visExport(type = "pdf", name = paste(group_id, desc, sep = "_"), label = "save") %>%
    addExport()
}
