#RC.clean.clean.gene.DESeqN.EB
dim(RC.clean.clean.gene.DESeqN.EB)

##############################

# 准备对EB数据进行remove embryo effect 操作

RC.clean.clean.gene.DESeqN.EB.big = 
  RC.clean.clean.gene.DESeqN.EB[rowMeans(RC.clean.clean.gene.DESeqN.EB) > 1,]

library(lme4) 

EB.embryo.info = factor(QC.clean.clean$embryo.number[QC.clean.clean$Stage=="E"])

fit.one.gene.EB <- function(x){
  b = lmer( x - mean(x) ~  0 + (1|factor(EB.embryo.info) ) )
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


#################################

# look for highly variable gene













