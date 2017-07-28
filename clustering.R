# 使用不同的研究得到的lineage marker,做聚类，看能否分成不同lineage

###########################################################

# early blastocyst 应该能分类，先看看EB

#####################

# 使用 Defining the three cell lineages of the human blastocyst by single-cell RNA-seq.
# table s3

Blakeley.tb3 = read.csv("Defining the three cell lineages of the human blastocyst by single-cell RNA-seq/TableS3.csv")

RC.clean.clean.gene.DESeqN.EB.Blakeley.tb3 = 
  RC.clean.clean.gene.DESeqN.EB[na.omit(match(Blakeley.tb3$Gene,rownames(RC.clean.clean.gene.DESeqN.EB))),]

palette.breaks <- seq(-3, 3, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)

heatmap( log2(as.matrix(RC.clean.clean.gene.DESeqN.EB.Blakeley.tb3)+1),trace = "none",density = "none",
         #Colv = as.dendrogram(complete.cluster),
         #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
         #ColSideColors =  c("grey","red")[factor(type)],
         col=color.palette,
         breaks = palette.breaks,
         scale = c("row"),
         dendrogram = "both")

## early blastocyst 分不了类






##########################################

# 试试 late blastocyst是不是可以分成三类，用堂复仇的数据

Tang.RPKM = read.csv("tangfuchou数据/nsmb.2660-S2.csv")

Tang.RPKM.lateblast = Tang.RPKM[,c(63:92)]


Tang.RPKM.lateblast.Blakeley.tb3 = 
  Tang.RPKM.lateblast[match(Blakeley.tb3$Gene,Tang.RPKM$Gene_ID),]


heatmap( log2(as.matrix(Tang.RPKM.lateblast.Blakeley.tb3)+1),trace = "none",density = "none",
         #Colv = as.dendrogram(complete.cluster),
         #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
         #ColSideColors =  c("grey","red")[factor(type)],
         col=color.palette,
         #breaks = palette.breaks,
         scale = c("row"),
         dendrogram = "both")






