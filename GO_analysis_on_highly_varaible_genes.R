##############################

# highly variable gene 和Lineage marker很多是重复的，highly enriched。


TE.1000 = read.table("TE.1000")
intersect(as.character(TE.1000$V1),sig.big.gene)


EPI.1000 = read.table("EPI.1000")
intersect(as.character(EPI.1000$V1),sig.big.gene)


PE.1000 = read.table("PE.1000")
intersect(as.character(PE.1000$V1),sig.big.gene)


ICM.TE.diff = read.table("ICM.TE.diff.gene.1000")
intersect(as.character(ICM.TE.diff$V1),sig.big.gene)

