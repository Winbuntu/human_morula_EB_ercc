library(gplots)

RC.clean.clean.gene.DESeqN.morula.sig.big.gene = 
  RC.clean.clean.gene.DESeqN.morula[match(sig.big.gene,rownames(RC.clean.clean.gene.DESeqN.morula)),]


palette.breaks <- seq(-1.5, 1.5, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)

png(file = "marula,highv gene, heatmap.jpeg",width = 1000,height = 800)
heatmap( log2(as.matrix(RC.clean.clean.gene.DESeqN.morula.sig.big.gene)+1),trace = "none",density = "none",
          #Colv = as.dendrogram(complete.cluster),
          #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
          #ColSideColors =  c("grey","red")[factor(type)],
          col=color.palette,
          breaks = palette.breaks,
          scale = c("none"),
          dendrogram = "both")
dev.off()

