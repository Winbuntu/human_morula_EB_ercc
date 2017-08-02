dim(RC.clean.gene.DeseqN)


gene_length.data = data.frame(name = rownames(RC),gene_length)


gene_length.data$gene_length[match(rownames(RC.morula.corrected),gene_length.data$name)]

# 计算RPKM，看三个cluster里面RPKM差异


# Rudi.RPM = sweep(RC,2,QC$Mapped_reads,FUN="/")*1000000
# Rudi.RPKM = sweep(Rudi.RPM,1,gene_length.data$gene_length,FUN="/")*1000
# 

RC.morula.corrected.quan.norm.RPKM = 
  sweep(RC.morula.corrected.quan.norm,1,gene_length.data$gene_length[match(rownames(RC.morula.corrected.quan.norm),gene_length.data$name)],
        FUN="/")*1000


####看分类？

RC.morula.corrected.quan.norm.RPKM.TE = 
  RC.morula.corrected.quan.norm.RPKM[ match(Fuchoutang.EPI$V1,rownames(RC.morula.corrected.quan.norm.RPKM)),]
group = cutree(complete.cluster,4)
res = t(apply(RC.morula.corrected.quan.norm.RPKM.TE,1, function(x){as.numeric(by((x), group, mean))}))[,1:3]
boxplot.matrix( log2(res+1) )

wilcox.test(log2(res+1)[,3],log2(res+1)[,2])

t.test(log2(res+1)[,3],log2(res+1)[,1])


######################


RC.morula.corrected.quan.norm.RPKM.TE = 
  RC.morula.corrected.quan.norm.RPKM[ match(M.marker.com$gene[M.marker.com$lineage=="TE"],rownames(RC.morula.corrected.quan.norm.RPKM)),]
group = cutree(complete.cluster,4)
res = t(apply(RC.morula.corrected.quan.norm.RPKM.TE,1, function(x){as.numeric(by((x), group, mean))}))[,1:3]
boxplot.matrix( log2(res+1) )

wilcox.test(log2(res+1)[,1],log2(res+1)[,3])

t.test(log2(res+1)[,2],log2(res+1)[,1])



###

# 随机基因

RC.morula.corrected.quan.norm.RPKM.TE = 
  RC.morula.corrected.quan.norm.RPKM[ 3000:3700,]
group = cutree(complete.cluster,4)
res = t(apply(RC.morula.corrected.quan.norm.RPKM.TE,1, function(x){as.numeric(by((x), group, mean))}))[,1:3]
boxplot.matrix( log2(res+1) )

wilcox.test(log2(res+1)[,2],log2(res+1)[,3])

t.test(log2(res+1)[,2],log2(res+1)[,1])