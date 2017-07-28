TE.1000 = read.table("TE.1000")
intersect(as.character(TE.1000$V1),sig.big.gene)


EPI.1000 = read.table("EPI.1000")
intersect(as.character(EPI.1000$V1),sig.big.gene)


PE.1000 = read.table("PE.1000")
intersect(as.character(PE.1000$V1),sig.big.gene)
