#' Parsing protein complex information from CORUM.
#' @description This function parse the dataset from CORUM (ver.3.0) to obtain protein membership to protein copmlexes.
#' @param organism Single character vector. Select from: "Human",Mouse","Pig","Bovine","Rat","Mammalia","Rabbit","Dog","Hamster"  "MINK"
#' @details The data is complied from CORUM The resulting output can be used for `ref` object in netvis() function.
#' @return list consists of:
#' \itemize{
#'   \item comp.tbl - data.frame. It is consisting from complexID(complex identifier given in CORUM database), group_id(complex name),varname(subunit).
#'   \item comp.list - list version of comp.tbl. Each list elements corresponds to a complex.
#'   \item comp.raw - data.frame. Raw data from CORUM database "allComplexes.txt" from http://mips.helmholtz-muenchen.de/corum/#download.
#' }
#' @references Giurgiu, M. et al. (2019) CORUM: the comprehensive resource of mammalian protein complexes-2019. Nucleic Acids Res, 47(D1), D559–D563.
#' @export
#' @examples
#' ref = parse.corum("Human")
#' names(ref)

parse.corum = function(organism = "Human"){
  data(corum.raw)
  cplx = corum.raw
  cplx = cplx[cplx$Organism == organism,]
  cplx_gene1 = strsplit(cplx$`subunits(Gene name)`,";|,| ")
  cplx_gene2 = strsplit(cplx$`subunits(Gene name syn)`,";|,| ")

  #cplx_gene = lapply(cplx_gene,function(x)unique(x[x!="None" & x!=""]))
  cplxMember1 = lapply(1:nrow(cplx),function(x)cbind(cplx$ComplexID[x],cplx$ComplexName[x],
                                                     unlist(cplx_gene1[x])))
  cplxMember2 = lapply(1:nrow(cplx),function(x)cbind(cplx$ComplexID[x],cplx$ComplexName[x],
                                                     unlist(cplx_gene2[x])))

  cplxMember1 = do.call(rbind,cplxMember1)
  cplxMember2 = do.call(rbind,cplxMember2)
  cplxMember = rbind(cplxMember1,cplxMember2)
  colnames(cplxMember) = c("group_id","desc","varname")
  cplxMember = cplxMember[cplxMember[,2]!="None"&cplxMember[,2]!="",]
  cplxMember = data.frame(cplxMember,check.names = FALSE,stringsAsFactors = FALSE)
  cplxMember = cplxMember[cplxMember$varname!="",]
  cplxMember = cplxMember[cplxMember$varname!="None",]
  cplxMember = cplxMember[order(cplxMember$group_id),]

  cplxMemberList = split(cplxMember,f = cplxMember$group_id)

  output = list(comp.tbl = cplxMember, comp.list = cplxMemberList, comp.raw = cplx)
  output
}
