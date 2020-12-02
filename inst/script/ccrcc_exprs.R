#require(readxl)
#require(data.table)


file_1 = system.file("extdata","CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv",package = "RoDiCE")

rawexpr =
  data.table::fread(file_1,
        header = F,skip=4,drop=1,sep = "\t",data.table = F)

sampleName =
  data.table::fread(file_1,
        header = F,nrows=1,drop=1,sep = "\t",data.table = F)

sampleName = unlist(sampleName)

geneName = data.table::fread(file_1,
                 header = F,skip=4,select=1,data.table = F)

geneName = unlist(geneName)

selectCol_1 = grep("Unshared",sampleName,invert = T)
selectCol_2 = grep("QC",sampleName,invert = T)
selectCol_3 = grep("^NCI",sampleName,invert = T)
selectCol_4 = grep("Log Ratio",sampleName)
selectCol = intersect(selectCol_1,intersect(selectCol_2,intersect(selectCol_3,selectCol_4)))

sampleName = sampleName[selectCol]
sampleName = gsub(" Log Ratio","",sampleName)

file_2 = system.file("extdata","S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.xlsx",package = "RoDiCE")
sampleSheet = read_excel::read_excel(file_2,
                         sheet=1,col_names = T)
sampleSheet = as.data.frame(sampleSheet)

exprList = list(expr = t(rawexpr[,selectCol]),
                geneName = geneName,
                sampleName = sampleName)


sampleSheet = sampleSheet[match(exprList$sampleName,sampleSheet$`Aliquot ID`),]
exprList$sampleGroup = factor(sampleSheet$Group)
exprList$pid = factor(sampleSheet$ParticipantID)

ccrcc = list()
ccrcc$tumor = exprList$expr[exprList$sampleGroup=="Tumor",]
ccrcc$normal = exprList$expr[exprList$sampleGroup=="Normal",]
colnames(ccrcc$tumor) = colnames(ccrcc$normal) = exprList$geneName

#library(pcaMethods)
pcaObj_tumor = pcaMethods::pca(ccrcc$tumor,nPcs=10)
pcaObj_normal = pcaMethods::pca(ccrcc$normal,nPcs=10)
ccrcc$tumor = pcaObj_tumor@completeObs
ccrcc$normal = pcaObj_normal@completeObs


file_3 = file.path(system.file("data",package = "RoDiCE"),"ccrcc.rda")
# file_3 = "/home/bmhi/projects/RoDiCE/data/ccrcc.rda"
save(ccrcc,file = file_3)

data(corum.hsp)
#6789,5273,6909,149,929
selectComp = corum.hsp[corum.hsp$group_id%in%c(6789,5273,6909,149,929),]
gene_names = unique(selectComp$varname)

ccrcc.subset = ccrcc
ccrcc.subset$tumor = ccrcc.subset$tumor[,colnames(ccrcc.subset$tumor)%in%gene_names]
ccrcc.subset$normal = ccrcc.subset$normal[,colnames(ccrcc.subset$normal)%in%gene_names]

file_4 = file.path(system.file("data",package = "RoDiCE"),"ccrcc.subset.rda")
# file_4 = "/home/bmhi/projects/RoDiCE/data/ccrcc.subset.rda"
save(ccrcc.subset,file = file_4)


selectComp = corum.hsp[corum.hsp$group_id%in%149,]
gene_names = unique(selectComp$varname)

ccrcc.pbaf = ccrcc
ccrcc.pbaf$tumor = ccrcc.pbaf$tumor[,colnames(ccrcc.pbaf$tumor)%in%gene_names]
ccrcc.pbaf$normal = ccrcc.pbaf$normal[,colnames(ccrcc.pbaf$normal)%in%gene_names]

file_5 = file.path(system.file("data",package = "RoDiCE"),"ccrcc.pbaf.rda")
# file_5 = "/home/bmhi/projects/RoDiCE/data/ccrcc.pbaf.rda"
save(ccrcc.pbaf,file = file_5)
