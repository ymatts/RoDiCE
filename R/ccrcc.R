#' Protein expression of clear renal cell carcinoma
#'
#' @description  Protein abudance dataset of clear renal cell carcinoma (ccRCC) published by CPTAC/TCGA (Clark et al., 2019). Data are available from the CPTAC data portal in the CPTAC Clear Cell Renal Cell Carcinoma (CCRCC) Discovery Study The data labeled "CPTAC_CompRef_CCRCC_Proteome_CDAP_Protein_Report.r1" were used.
#' @details A list of two data.frame of `tumor`(110 rows and 49 columns) and `normal`(110 rows and 49 columns). The rows and columns correspond to the sample and the protein, respectively
#' @format
#' @examples
#' data(ccrcc.subset) # Abudance matrix for 49 proteins as an example dataset
#' data(crccc.pbaf) # Abudance matrix for 10 proteins beloging to the PBAF complex as an example dataset
#'
#' # There are normal and tumor group samples in each element of the list.
#' dim(ccrcc.subset$tumor)
#' dim(ccrcc.subset$normal)
#' dim(ccrcc.pbaf$tumor)
#' dim(ccrcc.pbaf$normal)
#' @source \url{https://cptac-data-portal.georgetown.edu}
#' @references Clerk DJ et al.(2019) Integrated Proteogenomic Characterization of Clear Cell Renal Cell Carcinoma.Cell;179(4),964-983 e931.
#' @name ccrcc
NULL
#' @rdname ccrcc
"ccrcc.subset"
#' @rdname ccrcc
"ccrcc.pbaf"
