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


M.marker.com = read.table("maintained.markers.complete",header = T)
#intersect(as.character(Maintained.lineage.markers$V1),sig.big.gene)


Fuchoutang.TE = read.table("fuchoutang_TE")

Fuchoutang.EPI = read.table("fuchoutang.EPI")






TEAD4.target = read.table("TEAD4_target.gene.txt")
TEAD4.target$V1 = toupper(TEAD4.target$V1)

intersect(TEAD4.target$V1,sig.gene.after.correction)

intersect(TEAD4.target$V1,sig.big.gene)

####################################


phyper(length(intersect(as.character(PE.1000$V1),sig.big.gene)),
       length(sig.big.gene), dim(RC.clean.clean.gene.DESeqN.morula.big)[1]-length(sig.big.gene),
       length(as.character(PE.1000$V1)), lower.tail=TRUE);
#dim(RC.clean.clean.gene.DESeqN.morula.big)[1]
# (success-in-sample, success-in-bkgd, failure-in-bkgd, sample-size).

# 
# ############################################
# 
# # 画GO 气泡图
# 
# library(clusterProfiler)
# library(org.Hs.eg.db)
# 
# eg.morulaHigh.downreg.gene.DEseq = 
#   bitr(morulaHigh.downreg.gene.DEseq, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# 
# ggo.morulaHigh.downreg.gene.DEseq <- groupGO(gene     = eg.morulaHigh.downreg.gene.DEseq$ENTREZID,
#                OrgDb    = org.Hs.eg.db,
#                ont      = "CC",
#                level    = 3,
#                readable = TRUE)
# 
# head(ggo.morulaHigh.downreg.gene.DEseq)
# 
# ego.morulaHigh.downreg.gene.DEseq <- enrichGO(gene  = eg.morulaHigh.downreg.gene.DEseq$ENTREZID,
#                 #universe      = names(geneList),
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "CC",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)
# head(ego.morulaHigh.downreg.gene.DEseq)
# 
# dotplot(ego.morulaHigh.downreg.gene.DEseq)
# 
# ego2 <- enrichGO(gene         = eg.morulaHigh.downreg.gene.DEseq$ENTREZID,
#                  OrgDb         = org.Hs.eg.db,
#                  #keytype       = 'ENSEMBL',
#                  ont           = "BP",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 0.05
#                  #qvalueCutoff  = 0.05
#                  )
# dotplot(ego2)

#####################

# 自己画图

library(RColorBrewer)

morulaHigh.downreg.gene.DEseq.GO = read.csv("morulaLow.upreg.gene.ontology.csv")

ggplot(morulaHigh.downreg.gene.DEseq.GO[1:10,],
       aes(x=GO.biological.process.complete, y=upload_1..fold.Enrichment.,fill = upload_1..P.value. )) +
  geom_bar(stat="identity")+ theme_bw() + coord_flip() + scale_fill_gradient(low="red", high="blue")+
  theme(axis.text=element_text(size=16))

#print(p)
morulaHigh.downreg.gene.DEseq.GO = read.csv("morulaHigh.downreg.gene.DEseq.Geneontology.csv")

ggplot(morulaHigh.downreg.gene.DEseq.GO[1:10,],
       aes(x=GO.biological.process.complete, y=upload_1..fold.Enrichment.,fill = upload_1..P.value. )) +
  geom_bar(stat="identity")+ theme_bw() + coord_flip() + scale_fill_gradient(low="red", high="blue")+
  theme(axis.text=element_text(size=16))

