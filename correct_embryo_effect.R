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

#RPKM.clean.clean.norm.Morula.remove.EMB.eff = 2^RPKM.clean.clean.norm.Morula.log.remove.EMB.eff