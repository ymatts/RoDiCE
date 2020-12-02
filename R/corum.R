#' Protein membership to the protein complex
#'
#' @description  Protein membership information in protein complexes, based on the database CORUM(version 3.0).
#' @examples
#' data(corum.subset)
#' data(corum.hsp)
#' data(corum.hsp.pbaf)
#' data(corum.raw)
#' @details `corum.subset` is the example data used in vignette.
#' `corum.hsp` is the complex membership information of proteins restricted to "Human", and parse.corum("Human") provides the same.
#' `corum.hsp.pbaf` is the membership information of the PBAF complex only.
#' `corum.raw` is the raw data from CORUM.
#' \itemize{
#'   \item group_id. The protein complex identifier provided as `complexID` from CORUM.
#'   \item desc. The protein complex name provided as `complexName` from CORUM.
#'   \item varname. A group of proteins that belong to a protein complex.
#' }
#
#' @format
#' @source \url{http://mips.helmholtz-muenchen.de/corum/}
#' @references Giurgiu, M. et al. (2019) CORUM: the comprehensive resource of mammalian protein complexes-2019. Nucleic Acids Res, 47(D1), D559–D563.
#' @name corum
NULL
#' @rdname corum
"corum.raw"
#' @rdname corum
"corum.hsp"
#' @rdname corum
"corum.hsp.pbaf"
#' @rdname corum
"corum.subset"
