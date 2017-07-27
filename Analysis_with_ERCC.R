# R package loadin g start from here

library(data.table)
library(magrittr)
library(ggplot2)
library(ggthemes)
library(FactoMineR)
library(DESeq2)
library(corrplot)

###########
# 读取数据 QC information

QC = read.table("Human_EB_morula_QC_summery.tsv",header = F)
colnames(QC) = c("Sample","Total_reads","Mapped_reads","MT_reads")
# QC 表格的内容： ${TOTAL_READS}  ${MAPPED_READS} ${MT_READS}

QC$Fraction_mapped_reads = QC$Mapped_reads/QC$Total_reads
QC$MT_reads_fraction = QC$MT_reads/QC$Total_reads

QC$Stage = factor(substr(QC$Sample,1,1))

QC$gene.readcount = colSums(RC[-c(26365:26456),])

QC$gene.readcount.percent = QC$gene.readcount/QC$Mapped_reads

QC$embryo = sapply(strsplit(colnames(RC), split='-',fixed=TRUE),function(x) x[1])

QC$embryo.number = substr(QC$embryo,nchar(QC$embryo),nchar(QC$embryo))

############




###########
# 读取数据 Read count

RC = fread("Human_EB_Morula_READ_COUNT_table.txt")
gene_length = RC$Length

RC = as.data.frame(RC)
rownames(RC) = RC$Geneid
RC = RC[,-c(1:6)]
colnames(RC)
colnames(RC) = sapply(strsplit(colnames(RC), split='_',fixed=TRUE),function(x) x[1])

###### 
# 计算 RPM
Rudi.RPM = sweep(RC,2,QC$Mapped_reads,FUN="/")*1000000
QC$RPM.bigger.10.gene = apply(Rudi.RPM,2,function(x){  sum(x>10) } )
hist(QC$RPM.bigger.10.gene)

QC$detected.gene.number = apply(RC,2,function(x){  sum(x>0) } )


###########
# QC 图

# 
plot.mapped.reads <- 
  ggplot(QC, aes(x = Total_reads, y = Fraction_mapped_reads , shape = (Stage), colour = (Stage) )) + 
  geom_point(size = 4) + theme_base()
print(plot.mapped.reads)

# 
plot.F.MT.reads <- 
  ggplot(QC, aes(x = Total_reads, y = MT_reads_fraction , shape = (Stage), colour = (Stage) )) + 
  geom_point(size = 4) + theme_base()
print(plot.F.MT.reads)

#

plot.gene.10RPM.reads <- 
  ggplot(QC, aes(x = Total_reads, y = RPM.bigger.10.gene , shape = (Stage), colour = (Stage) )) + 
  geom_point(size = 4) + theme_base()
print(plot.gene.10RPM.reads)

#
plot.gene.reads <- 
  ggplot(QC, aes(x = Total_reads, y = gene.readcount.percent , shape = (Stage), colour = (Stage) )) + 
  geom_point(size = 4) + theme_base()
print(plot.gene.reads)


#

pca.res = PCA(QC[,c(5,6,7,10)],graph = F)


PCs <- data.frame(PC1 = as.numeric(pca.res$ind$coord[,1]) ,
                  PC2 = as.numeric(pca.res$ind$coord[,2]))

plot(PCs$PC1,PCs$PC2)
abline(v=-3,col="red")

####################################

PCA.bad.cell = PCs$PC1 < -8
QC[PCA.bad.cell,]

#######

Gene.reads.filter = QC$gene.readcount < 500000 # 暂时选择500,000  之后可以提高到1M

sample.filter = (Gene.reads.filter | PCA.bad.cell)


######################

# 这里删除不好的sample
RC.clean = RC[,!sample.filter]
QC.clean = QC[!sample.filter,]


########################################################

RC.clean.ERCC = RC.clean[c(26365:26456),]
RC.clean.gene = RC.clean[-c(26365:26456),]


# 开始使用Deseq2 normalize sample

cts.ERCC = RC.clean.ERCC

coldata = data.frame( stage = QC.clean[,c(8)])
rownames(coldata) = colnames(cts.ERCC)

#coldata$stage = factor(c("E",rep("M",169))) # 测试发现，sizefactor并不随着group信息而变化

dds.ERCC <- DESeqDataSetFromMatrix(countData = cts.ERCC,
                              colData = coldata,
                              design = ~ stage)

dds.ERCC <- estimateSizeFactors(dds.ERCC)

plot(colSums(RC.clean.ERCC) ,colData(dds.ERCC)$sizeFactor,xlim= c(0,7000000))
plot(colSums(counts(dds.ERCC, normalized = TRUE)) ,colData(dds.ERCC)$sizeFactor,xlim= c(0,7000000))


RC.clean.ERCC.DeseqN = (counts(dds.ERCC, normalized = TRUE))

####################

# 使用Deseq2 normalize gene reads

cts.gene = RC.clean.gene

dds.Gene <- DESeqDataSetFromMatrix(countData = cts.gene,
                                   colData = coldata,
                                   design = ~ stage) 

dds.Gene <- estimateSizeFactors(dds.Gene)

plot(colSums(RC.clean.gene) ,colData(dds.Gene)$sizeFactor,xlim= c(0,7000000))

plot(colSums(counts(dds.Gene, normalized = TRUE)) ,colData(dds.Gene)$sizeFactor,xlim= c(0,7000000))

RC.clean.gene.DeseqN =  counts(dds.Gene, normalized = TRUE)

#########################


GetPCA.Norm.data <- function(RPM.data){
  #RPM.data = RPM.clean
  require(FactoMineR)
  gene_filter  = apply(   RPM.data  ,1, function(x) {sum(x>10) >= ( dim(RPM.data)[2]  /10)   } )
  RPM.data.filtered = RPM.data[ gene_filter ,]
  
  require(preprocessCore)
  RPM.data.filtered.norm = normalize.quantiles.robust(  (as.matrix(RPM.data.filtered))   ,use.median=TRUE,use.log2=FALSE)
  colnames(RPM.data.filtered.norm) = colnames(RPM.data.filtered)
  rownames(RPM.data.filtered.norm) = rownames(RPM.data.filtered)
  
  
###########
  
  pca.res = PCA(t(  log2(  RPM.data.filtered.norm +0.5)  ),graph = F)
  
  
  PCs <- data.frame(PC1 = as.numeric(pca.res$ind$coord[,1]) ,
                    PC2 = as.numeric(pca.res$ind$coord[,2]))
  
  
  return(list(PCs, RPM.data.filtered.norm))
}


#### 

# 使用Gene read 画PCA
total.pca.plot.Gene = ggplot(GetPCA.Norm.data(RC.clean.gene.DeseqN)[[1]], 
                        aes(PC1,PC2,shape = factor(QC.clean$Stage),
                            colour = factor(QC.clean$embryo.number),
                            size = (colData(dds.Gene)$sizeFactor)
                            )) + 
  geom_point() + theme_base() 
print(total.pca.plot.Gene)

#

total.pca.plot.Gene = ggplot(GetPCA.Norm.data(RC.clean.gene.DeseqN)[[1]], 
                             aes(PC1,PC2,shape = factor(QC.clean$Stage),
                                 colour = factor(QC.clean$embryo.number),
                                 size = (QC.clean$RPM.bigger.10.gene)
                             )) + 
  geom_point() + theme_base() 
print(total.pca.plot.Gene)


# 使用ERCC 画PCA

total.pca.plot.ERCC = ggplot(GetPCA.Norm.data(RC.clean.ERCC.DeseqN)[[1]], 
                        aes(PC1,PC2,shape = factor(QC.clean$Stage),
                            colour = factor(QC.clean$embryo.number),
                            size = (colData(dds.Gene)$sizeFactor)
                        )) + 
  geom_point() + theme_base() 
print(total.pca.plot.ERCC)

#############

# 检测到的基因数目，在两类细胞中是不是有差异？

# 使用bar plot，进行QC
a = GetPCA.Norm.data(RC.clean.gene.DeseqN)[[1]]

a$PC1 < (-50)

boxplot( list(QC.clean$RPM.bigger.10.gene[a$PC1 < (-50)],
              QC.clean$RPM.bigger.10.gene[a$PC1 > (-50)]
              ))

boxplot( list(QC.clean$detected.gene.number[a$PC1 < (-50)],
              QC.clean$detected.gene.number[a$PC1 > (-50)]
))


boxplot( list(QC.clean$gene.readcount.percent[a$PC1 < (-50)],
              QC.clean$gene.readcount.percent[a$PC1 > (-50)]
))


##########################################

# 之后去除左边subgroup的细胞

RC.clean.clean.gene.DESeqN = RC.clean.gene.DeseqN[,a$PC1 > (-50)]
RC.clean.clean.ERCC.DESeqN = RC.clean.ERCC.DeseqN[,a$PC1 > (-50)]

QC.clean.clean = QC.clean[a$PC1 > (-50),]


# 之后再PCA

ggplot(GetPCA.Norm.data(RC.clean.clean.gene.DESeqN)[[1]], 
       aes(PC1,PC2,shape = factor(QC.clean.clean$Stage),
           colour = factor(QC.clean.clean$embryo.number))) + 
  geom_point(size=3) + theme_base() 


##############

GetPCA.Norm.data.noquantile.norm <- function(RPM.data){
  #RPM.data = RPM.clean
  require(FactoMineR)
  gene_filter  = apply(   RPM.data  ,1, function(x) {sum(x>10) >= ( dim(RPM.data)[2]  /10)   } )
  RPM.data.filtered = RPM.data[ gene_filter ,]
  
  #require(preprocessCore)
  #RPM.data.filtered.norm = normalize.quantiles.robust(  (as.matrix(RPM.data.filtered))   ,use.median=TRUE,use.log2=FALSE)
  #colnames(RPM.data.filtered.norm) = colnames(RPM.data.filtered)
  #rownames(RPM.data.filtered.norm) = rownames(RPM.data.filtered)
  
  
  ###########
  
  pca.res = PCA(t(  log2(  RPM.data.filtered +1)  ),graph = F)
  
  
  PCs <- data.frame(PC1 = as.numeric(pca.res$ind$coord[,1]) ,
                    PC2 = as.numeric(pca.res$ind$coord[,2]))
  
  
  return(list(PCs, RPM.data.filtered))
}

#####

ggplot(GetPCA.Norm.data.noquantile.norm(RC.clean.clean.gene.DESeqN)[[1]], 
       aes(PC1,PC2,shape = factor(QC.clean.clean$Stage),
           colour = factor(QC.clean.clean$embryo.number))) + 
  geom_point(size=3) + theme_base() 





ggplot(GetPCA.Norm.data(RC.clean.clean.ERCC.DESeqN)[[1]], 
       aes(PC1,PC2,shape = factor(QC.clean.clean$Stage),
           colour = factor(QC.clean.clean$embryo.number))) + 
  geom_point(size=3) + theme_base() 


#########

# 对样品的ERCC作图，发现correlation都很高

corrplot(cor((RC.clean.clean.ERCC.DESeqN) ),method="color",order="hclust")


########

# 计算胚胎内／胚胎间的correlation

# morula stage 
RC.clean.clean.gene.DESeqN.morula = RC.clean.clean.gene.DESeqN[,QC.clean.clean$Stage=="M"]

coor.matrix.RC.clean.clean.gene.DESeqN.morula = 
  cor(log2(  RC.clean.clean.gene.DESeqN.morula  +1),method = "spearman")


sample.annotation.morula = QC.clean.clean$embryo.number[QC.clean.clean$Stage == "M"]

embryos.morula = unique(sample.annotation.morula)

intra.marula.corr = NULL

for(i in embryos.morula){
  position = print(sample.annotation.morula == i)
  one.embryo.corr.matrix = coor.matrix.RC.clean.clean.gene.DESeqN.morula[position,position]
  intra.cor = one.embryo.corr.matrix[lower.tri(one.embryo.corr.matrix, diag = FALSE)]
  intra.marula.corr = c(intra.marula.corr,intra.cor)
}


corr.ma.clean.clean.norm.morula.inter = coor.matrix.RC.clean.clean.gene.DESeqN.morula

for(i in embryos.morula){
  position = print(sample.annotation.morula == i)
  corr.ma.clean.clean.norm.morula.inter[position,position] = NA

}

morula.cor.data = data.frame(cor = c(as.numeric(corr.ma.clean.clean.norm.morula.inter),as.numeric(intra.marula.corr)),
                             label = c(rep("Inter",length(as.numeric(corr.ma.clean.clean.norm.morula.inter))) ,
                                       rep("Intra",length(as.numeric(intra.marula.corr)))  ) )

ggplot(morula.cor.data, aes(x = label, y = cor,fill = label)) +
  geom_boxplot()+theme_classic()

t.test(as.numeric(corr.ma.clean.clean.norm.morula.inter),as.numeric(intra.marula.corr))
# p-value < 2.2e-16

###########################
# EB stage 
RC.clean.clean.gene.DESeqN.EB = RC.clean.clean.gene.DESeqN[,QC.clean.clean$Stage=="E"]

corr.ma.clean.clean.norm.EB = cor(log2(  RC.clean.clean.gene.DESeqN.EB  +1),method = "spearman")


sample.annotation.EB = QC.clean.clean$embryo.number[QC.clean.clean$Stage == "E"]
#lower.tri(x, diag = FALSE)
embryos.EB = unique(sample.annotation.EB)

#lower.tri(corr.ma.clean.clean.norm.morula, diag = FALSE)

intra.EB.corr = NULL

for(i in embryos.EB){
  position = print(sample.annotation.EB == i)
  one.embryo.corr.matrix = corr.ma.clean.clean.norm.EB[position,position]
  intra.cor = one.embryo.corr.matrix[lower.tri(one.embryo.corr.matrix, diag = FALSE)]
  intra.EB.corr = c(intra.EB.corr,intra.cor)
}


corr.ma.clean.clean.norm.EB.inter = corr.ma.clean.clean.norm.EB

for(i in embryos.EB){
  position = print(sample.annotation.EB == i)
  corr.ma.clean.clean.norm.EB.inter[position,position] = NA
  #intra.cor = one.embryo.corr.matrix[lower.tri(one.embryo.corr.matrix, diag = FALSE)]
  #intra.marula.corr = c(intra.marula.corr,intra.cor)
}

EB.cor.data = data.frame(cor = c(as.numeric(corr.ma.clean.clean.norm.EB.inter),as.numeric(intra.EB.corr)),
                         label = c(rep("Inter",length(as.numeric(corr.ma.clean.clean.norm.EB.inter))) ,
                                   rep("Intra",length(as.numeric(intra.EB.corr)))  ) )

ggplot(EB.cor.data, aes(x = label, y = cor,fill = label)) +
  geom_boxplot()+theme_classic()


t.test(as.numeric(corr.ma.clean.clean.norm.EB.inter),as.numeric(intra.EB.corr))


##########

cor.data = rbind(morula.cor.data,EB.cor.data)
cor.data$stage = c(rep("M",dim(morula.cor.data)[1]),rep("E",dim(EB.cor.data)[1]))

inter.intra.coorlation.box.plot = ggplot(cor.data, aes(x = factor(stage,levels = c("M","E")), y = cor,
                                                       fill = factor(  as.character(label)  ,levels = c("Intra","Inter"))  )) +
  geom_boxplot()+theme_classic() + ylim(c(0.4,1))

print(inter.intra.coorlation.box.plot)



##############

#dds.ERCC
#coldata

dds.ERCC.50reads <- dds.ERCC[ rowMeans(counts(dds.ERCC)) > 50, ]
dds.ERCC.50reads <- DESeq(dds.ERCC.50reads)


res.ERCC.50reads <- results(dds.ERCC.50reads)
2^res.ERCC.50reads$log2FoldChange
sum(res.ERCC.50reads$padj < 0.05)
#


# ERCC 中没有差异表达的基因


####
# 看一下morula 和 blastocyst间的差异表达的基因

dds.Gene.50reads = dds.Gene[ rowMeans(counts(dds.Gene)) > 20, ]
dds.Gene.50reads <- DESeq(dds.Gene.50reads)
res.Gene.50reads <- results(dds.Gene.50reads)

## morula 中高表达,随时间下调的基因

sum((res.Gene.50reads$log2FoldChange > log2(2)) & ( res.Gene.50reads$padj < 0.05))

morulaHigh.downreg.gene.DEseq = 
  rownames(res.Gene.50reads)[((res.Gene.50reads$log2FoldChange > log2(2)) & ( res.Gene.50reads$padj < 0.05))]

morulaLow.upreg.gene.DEseq = 
  rownames(res.Gene.50reads)[((res.Gene.50reads$log2FoldChange < log2(1/2)) & ( res.Gene.50reads$padj < 0.05))]

write.table(morulaHigh.downreg.gene.DEseq,file = "morulaHigh.downreg.gene.DEseq.txt",quote = F,row.names = F)

write.table(morulaLow.upreg.gene.DEseq,file = "morulaLow.upreg.gene.DEseq.txt",quote = F,row.names = F)







