#' Protein membership to the protein complex
#'
#' @description  Protein membership information in protein complexes, based on the database CORUM(version 3.0).
#' @examples
#' data(corum.subset)
#' data(corum.hsp)
#' data(corum.hsp.pbaf)
#' @details `corum.subset` is the example data used in vignette.
#' `corum.hsp` is the complex membership information of proteins restricted to "Human", and parse.corum("Human") provides the same.
#' `corum.hsp.pbaf` is the membership information of the PBAF complex only.
#' \itemize{
#'   \item group_id. The protein complex identifier provided as `complexID` from CORUM.
#'   \item desc. The protein complex name provided as `complexName` from CORUM.
#'   \item varname. A group of proteins that belong to a protein complex.
#' }
#
#' @format
#' @source \url{https://cptac-data-portal.georgetown.edu}
#' @references Clerk DJ et al.(2019) Integrated Proteogenomic Characterization of Clear Cell Renal Cell Carcinoma.Cell;179(4),964-983 e931.
#' @name corum
NULL
#' @rdname corum
"corum.hsp"
#' @rdname corum
"corum.hsp.pbaf"
#' @rdname corum
"corum.subset"
