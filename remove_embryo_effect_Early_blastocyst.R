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
# 用之前记得从新运行，覆盖掉旧数据

RC.clean.clean.ERCC.DESeqN.EB = 
  RC.clean.clean.ERCC.DESeqN[,match(colnames(RC.clean.clean.gene.DESeqN.EB),
                                    colnames(RC.clean.clean.ERCC.DESeqN))]

meansHeLa <- rowMeans( RC.clean.clean.ERCC.DESeqN.EB )
varsHeLa <- rowVars( RC.clean.clean.ERCC.DESeqN.EB )
cv2HeLa <- varsHeLa / meansHeLa^2


minMeanForFit <- unname( quantile( meansHeLa[ which( cv2HeLa > .3 ) ], .95 ) )
minMeanForFit

useForFit <- (meansHeLa >= minMeanForFit)

fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansHeLa[useForFit] ),
                   cv2HeLa[useForFit] )
fit$coefficients

### 这里有错
sfHeLa = 
  colData(dds.ERCC)$sizeFactor[
    match(colnames(RC.clean.clean.gene.DESeqN.EB),names(colData(dds.ERCC)$sizeFactor))]

xi <- mean( 1 / sfHeLa )


a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"] - xi )


c( a0, a1 )

####

plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .0005, 100 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)")
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-3:2), c("0.001", "0.01", "0.1", "1", "10","100" ), las=2 )
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# Add the data points
points( meansHeLa, cv2HeLa, pch=20, cex=1, col="blue" )
# Plot the fitted curve
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
# Plot quantile lines around the fit


######
df <- ncol(RC.clean.clean.gene.DESeqN.morula) - 1  # 
lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .975, df ) / df,
       col="#FF000080", lwd=2, lty="dashed" )
lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .025, df ) / df,
       col="#FF000080", lwd=2, lty="dashed" )









