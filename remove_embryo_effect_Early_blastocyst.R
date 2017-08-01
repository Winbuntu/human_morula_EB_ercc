#RC.clean.clean.gene.DESeqN.EB
dim(RC.clean.clean.gene.DESeqN.EB)

##############################

# 准备对EB数据进行remove embryo effect 操作

RC.clean.clean.gene.DESeqN.EB.big = 
  RC.clean.clean.gene.DESeqN.EB[rowMeans(RC.clean.clean.gene.DESeqN.EB) > 1,]


###########################################

RC.clean.clean.gene.DESeqN.EB.maintained.marker = 
  RC.clean.clean.gene.DESeqN.EB[na.omit(match(Maintained.lineage.markers$V1,rownames(RC.clean.clean.gene.DESeqN.EB))),]

heatmap( log2(as.matrix(RC.clean.clean.gene.DESeqN.EB.maintained.marker)+1),trace = "none",density = "none",
         #Colv = as.dendrogram(complete.cluster),
         #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
         #ColSideColors =  c("grey","red")[factor(type)],
         col=color.palette,
         breaks = palette.breaks,
         scale = c("row"),
         dendrogram = "both")




library(lme4) 

EB.embryo.info = factor(QC.clean.clean$embryo.number[QC.clean.clean$Stage=="E"])

fit.one.gene.EB <- function(x){
  b = lmer( x - mean(x) ~   (1|factor(EB.embryo.info) ) )
  return(x - fitted(b))
}


RC.clean.clean.gene.DESeqN.EB.big.log = log2(RC.clean.clean.gene.DESeqN.EB.big+1)

RC.clean.clean.gene.DESeqN.EB.big.log.removed.Embryo.effect = 
  apply(RC.clean.clean.gene.DESeqN.EB.big.log,1,fit.one.gene.EB)

RC.clean.clean.gene.DESeqN.EB.big.log.removed.Embryo.effect = t(RC.clean.clean.gene.DESeqN.EB.big.log.removed.Embryo.effect)

RC.clean.clean.gene.DESeqN.EB.big.NOT.log.removed.Embryo.effect = 
  2^RC.clean.clean.gene.DESeqN.EB.big.log.removed.Embryo.effect-1


##################################

# 画PCA plot看差异

ggplot(GetPCA.Norm.data.noquantile.norm(RC.clean.clean.gene.DESeqN.EB.big.NOT.log.removed.Embryo.effect)[[1]], 
       aes(PC1,PC2,shape = EB.embryo.info,
           colour = EB.embryo.info)) + 
  geom_point(size=3) + theme_base() 

##

ggplot(GetPCA.Norm.data.noquantile.norm(RC.clean.clean.gene.DESeqN.EB)[[1]], 
       aes(PC1,PC2,shape = EB.embryo.info,
           colour = EB.embryo.info)) + 
  geom_point(size=3) + theme_base() 


##########################################################

# 使用maintained marker无法聚类

RC.clean.clean.gene.DESeqN.EB.big.NOT.log.removed.Embryo.effect.maintained.marker = 
  RC.clean.clean.gene.DESeqN.EB.big.NOT.log.removed.Embryo.effect[na.omit(match(Maintained.lineage.markers$V1,
                                                                                rownames(RC.clean.clean.gene.DESeqN.EB.big.NOT.log.removed.Embryo.effect))),]

heatmap( log2(as.matrix(RC.clean.clean.gene.DESeqN.EB.big.NOT.log.removed.Embryo.effect.maintained.marker)+1),trace = "none",density = "none",
         #Colv = as.dendrogram(complete.cluster),
         #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
         #ColSideColors =  c("grey","red")[factor(type)],
         col=color.palette,
         breaks = palette.breaks,
         scale = c("row"),
         dendrogram = "both")



#################################

# look for highly variable gene








