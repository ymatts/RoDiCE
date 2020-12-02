#library(data.table)
infile_1 = system.file("extdata","allComplexes.txt",package = "RoDiCE")
corum.raw = fread(infile_1,data.table=FALSE)

outfile_0 = file.path("corum.raw",system.file("data",package = "RoDiCE"))
# outfile_0 = "/home/bmhi/projects/RoDiCE/data/corum.raw.rda"
save(corum.raw,file=outfile_0)

outfile_1 = file.path("corum.hsp",system.file("data",package = "RoDiCE"))
# outfile_1 = "/home/bmhi/projects/RoDiCE/data/corum.hsp.rda"
corum_obj = parse.corum("Human")
corum.hsp = corum_obj$comp.tbl
save(corum.hsp,file=outfile_1)

outfile_2 = file.path("corum.hsp.pbaf",system.file("data",package = "RoDiCE"))
# outfile_2 = "/home/bmhi/projects/RoDiCE/data/corum.hsp.pbaf.rda"
corum.hsp.pbaf = corum.hsp[corum.hsp$group_id==149,]
save(corum.hsp.pbaf,file=outfile_2)


outfile_3 = file.path("corum.subset",system.file("data",package = "RoDiCE"))
# outfile_3 = "/home/bmhi/projects/RoDiCE/data/corum.subset.rda"
corum_obj = parse.corum("Human")
corum.hsp = corum_obj$comp.tbl
corum.subset = corum.hsp[corum.hsp$group_id %in%c(6789,5273,6909,149,929),]
save(corum.subset,file=outfile_3)
