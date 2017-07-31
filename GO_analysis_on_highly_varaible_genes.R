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

Maintained.lineage.markers = read.table("cell.matintained.3.lineage.markers")
intersect(as.character(Maintained.lineage.markers$V1),sig.big.gene)


TEAD4.target = read.table("TEAD4_target.gene.txt")
TEAD4.target$V1 = toupper(TEAD4.target$V1)

intersect(TEAD4.target$V1,sig.gene.after.correction)

intersect(TEAD4.target$V1,sig.big.gene)

