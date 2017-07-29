###############################

# 使用linear mixed model 去除胚胎effect 
# 先选择表达大于1的基因
RC.clean.clean.gene.DESeqN.morula.big = RC.clean.clean.gene.DESeqN.morula[rowMeans(RC.clean.clean.gene.DESeqN.morula) > 1,]

library(lme4) 

morula.embryo.info = factor(QC.clean.clean$embryo.number[QC.clean.clean$Stage=="M"])

fit.one.gene <- function(x){
  b = lmer( x - mean(x) ~  0 + (1|factor(morula.embryo.info) ) )
  return(x - fitted(b))
}

RC.clean.clean.gene.DESeqN.morula.big.log = log2(RC.clean.clean.gene.DESeqN.morula.big+1)

RC.clean.clean.gene.DESeqN.morula.big.log.removed.Embryo.effect = apply(RC.clean.clean.gene.DESeqN.morula.big.log,1,fit.one.gene)

RC.clean.clean.gene.DESeqN.morula.big.log.removed.Embryo.effect = t(RC.clean.clean.gene.DESeqN.morula.big.log.removed.Embryo.effect)

RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect = 
  2^RC.clean.clean.gene.DESeqN.morula.big.log.removed.Embryo.effect-1

#RPKM.clean.clean.norm.Morula.remove.EMB.eff = 2^RPKM.clean.clean.norm.Morula.log.remove.EMB.eff

ggplot(GetPCA.Norm.data.noquantile.norm(RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect)[[1]], 
       aes(PC1,PC2,shape = morula.embryo.info,
           colour = morula.embryo.info)) + 
  geom_point(size=3) + theme_base() 


ggplot(GetPCA.Norm.data.noquantile.norm(RC.clean.clean.gene.DESeqN.morula)[[1]], 
       aes(PC1,PC2,shape = morula.embryo.info,
           colour = morula.embryo.info)) + 
  geom_point(size=3) + theme_base() 

##############################################

RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect

require(FactoMineR)

RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA = 
  RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect[ apply(   RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect  ,1, 
                                              function(x) {sum(x>10) >= ( dim(RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect)[2]  /10)   } ) ,]

res.PCA.RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA = 
  PCA(t(  log2(  RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA +1)  ),graph = F)

PCs.corrected <- data.frame(PC1 = as.numeric(res.PCA.RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA$ind$coord[,1]) ,
                  PC2 = as.numeric(res.PCA.RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA$ind$coord[,2]))

# 这是经过修正的PCA

ggplot(PCs.corrected, 
       aes(PC1,PC2,shape = morula.embryo.info,
           colour = morula.embryo.info)) + 
  geom_point(size=3) + theme_base() 


pc.loading.corrected = 
  res.PCA.RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA$var$coord[,1]^2 +  
  res.PCA.RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA$var$coord[,2]^2

pc.loading.corrected.ordered = pc.loading.corrected[order(pc.loading.corrected,decreasing = T)]


res.PCA.morula.corrected.ordered = 
  res.PCA.RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA$var$coord[order(pc.loading.corrected,decreasing = T),]

length(Blakeley.tb3$Gene)

intersect(names(pc.loading.corrected.ordered)[1:1000],Blakeley.tb3$Gene)

intersect(rownames(res.PCA.morula.corrected.ordered)[1:1000],Blakeley.tb3$Gene)

#######
# 发现有更多的基因重合了

plot(res.PCA.morula.corrected.ordered[1:1000,1],
     res.PCA.morula.corrected.ordered[1:1000,2],pch="."  )

TE.pos = na.omit(match(TE.1000$V1,rownames(res.PCA.morula.corrected.ordered)))

TE.pos[TE.pos>1000] <- NA

points(res.PCA.morula.corrected.ordered[TE.pos,1],
       res.PCA.morula.corrected.ordered[TE.pos,2],col="red")


################

EPI.pos = na.omit(match(EPI.1000$V1,rownames(res.PCA.morula.corrected.ordered)))

EPI.pos[EPI.pos>1000] <- NA

points(res.PCA.morula.corrected.ordered[EPI.pos,1],
       res.PCA.morula.corrected.ordered[EPI.pos,2],col="blue")

##################

PE.pos = na.omit(match(PE.1000$V1,rownames(res.PCA.morula.corrected.ordered)))

PE.pos[PE.pos>1000] <- NA

points(res.PCA.morula.corrected.ordered[PE.pos,1],
       res.PCA.morula.corrected.ordered[PE.pos,2],col="green")



##################

# 看top PC 的基因，之间的correlation 是怎样的

RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA.top1000 = 
  RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA[match(names(pc.loading.corrected.ordered)[1:1000],rownames(RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA)),]

morula.corrected.top1000.cor = cor( t(
  
  log2(RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect.filtered.for.PCA.top1000+1)
  
  
  ),method = "spearman")

#corrplot( morula.corrected.top1000.cor, method="color",tl.pos = "r")
palette.breaks <- seq(-1, 1, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)

marker.color = ifelse( rownames(morula.corrected.top1000.cor) %in%  TE.1000$V1, "red","white" )
marker.color[rownames(morula.corrected.top1000.cor) %in%  EPI.1000$V1] = "blue"
#marker.color[rownames(morula.corrected.top1000.cor) %in%  $V1] = "blue"

png(file = "marula,highv gene, correlation matrix heatmap.jpeg",width = 1000,height = 800)
heatmap( morula.corrected.top1000.cor ,trace = "none",density = "none",
         #Colv = as.dendrogram(complete.cluster),
         #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
         ColSideColors =  marker.color,
         col=color.palette,
         breaks = palette.breaks,
         scale = c("none"),
         dendrogram = "both")
dev.off()














