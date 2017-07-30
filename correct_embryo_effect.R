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



#####################

# 试着找出top variable gene，在embryo effect corrected的data里

# 用之前要重新运行，以免被覆盖


# RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect

RC.morula.corrected = RC.clean.clean.gene.DESeqN.morula.big.NOT.log.removed.Embryo.effect


##### 

# technique noise


RC.clean.clean.ERCC.DESeqN.Morula = 
  RC.clean.clean.ERCC.DESeqN[,match(colnames(RC.clean.clean.gene.DESeqN.morula),
                                    colnames(RC.clean.clean.ERCC.DESeqN))]

meansHeLa <- rowMeans( RC.clean.clean.ERCC.DESeqN.Morula )
varsHeLa <- rowVars( RC.clean.clean.ERCC.DESeqN.Morula )
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
    match(colnames(RC.clean.clean.gene.DESeqN.morula),names(colData(dds.ERCC)$sizeFactor))]

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




##############################################
# perform actural test

meansAt <- rowMeans( RC.morula.corrected )
varsAt <- rowVars( RC.morula.corrected )
cv2At <- varsAt / meansAt^2

sfAt = 
  colData(dds.Gene)$sizeFactor[
    match(colnames(RC.morula.corrected),names(colData(dds.Gene)$sizeFactor))]

psia1theta <- mean( 1 / sfAt ) + a1 * mean( sfHeLa / sfAt )

minBiolDisp <- .5^2

m <- ncol(RC.morula.corrected)
cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
testDenom <- ( meansAt * psia1theta + meansAt^2 * cv2th ) / ( 1 + cv2th/m )
p <- 1 - pchisq( varsAt * (m-1) / testDenom, m-1 )

padj <- p.adjust( p, "BH" )
sig <- padj < .1
sig[is.na(sig)] <- FALSE
table( sig )


plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)")
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10","100" ), las=2 )
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# Plot the plant genes, use a different color if they are highly variable
points( meansAt, cv2At, pch=20, cex=.2,
        col = ifelse( padj < .1, "#C0007090", colAt ) )
# Add the technical noise fit, as before
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
# Add a curve showing the expectation for the chosen biological CV^2 thershold
lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3 )

abline(v = 1)

sig.gene.after.correction = names(sig)[sig]
write.table(sig.gene.after.correction,file = "morula.high.variable.AFTERCorrectEmbryoEffect.txt",quote = F,row.names = F)      


intersect(sig.gene.after.correction,Blakeley.tb3$Gene)

intersect(sig.gene.after.correction,sig.big.gene)

###########

morula.gene.distance = cv2At - meansAt

morula.gene.distance.big = morula.gene.distance[meansAt>10]

morula.gene.distance.order = morula.gene.distance.big[order(morula.gene.distance.big,decreasing = T)]

match(names(morula.gene.distance.order[1:1000]),names(meansAt))

points( meansAt[match(names(morula.gene.distance.order[1:1000]),names(meansAt))], cv2At[match(names(morula.gene.distance.order[1:1000]),names(meansAt))], pch=20, cex=.2,
        col = "black" )

#morula.gene.distance.order[1:20]

padj.order = padj[order(padj,decreasing = F)]

RC.morula.corrected.top.variable = RC.morula.corrected[match(names(padj.order)[1:1000],rownames(RC.morula.corrected)),]

RC.morula.corrected.big = RC.morula.corrected.top.variable[rowMeans(RC.morula.corrected.top.variable) > 10,]



intersect(rownames(RC.morula.corrected.big),Blakeley.tb3$Gene)
intersect(rownames(RC.morula.corrected.big),TE.1000$V1)

palette.breaks <- seq(-2, 2, 0.1)
color.palette = colorRampPalette(c("dodgerblue4","dodgerblue1","white","firebrick1","firebrick3"), 
                                 space="Lab")(length(palette.breaks) - 1)


png(file = "marula,highv gene after correcting for embryo effect, heatmap.jpeg",width = 1000,height = 800)
heatmap( log2(as.matrix(RC.morula.corrected.big)+1),trace = "none",density = "none",
         #Colv = as.dendrogram(complete.cluster),
         #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
         #ColSideColors =  c("grey","red")[factor(type)],
         col=color.palette,
         breaks = palette.breaks,
         scale = c("row"),
         dendrogram = "both")
dev.off()

############################


plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)")
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10","100" ), las=2 )
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# Plot the plant genes, use a different color if they are highly variable
points( meansAt, cv2At, pch=20, cex=.2,
        col = ifelse( padj < .1, "#C0007090", colAt ) )
# Add the technical noise fit, as before
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
# Add a curve showing the expectation for the chosen biological CV^2 thershold
lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3 )

abline(v = 1)

points( meansAt[match(rownames(RC.morula.corrected.top.variable),names(meansAt))],
        cv2At[match(rownames(RC.morula.corrected.top.variable),names(meansAt))], pch=20, cex=.2,
        col = "black" )


###############################

#RC.morula.corrected.top.variable

png(file = "marula,highv gene after correcting for embryo effect, heatmap, including genes >1 count.jpeg",width = 1000,height = 800)
heatmap( log2(as.matrix(RC.morula.corrected.top.variable)+1),trace = "none",density = "none",
         #Colv = as.dendrogram(complete.cluster),
         #ColSideColors = c("grey","red")[ factor(type[complete.cluster$order] ) ] ,
         #ColSideColors =  c("grey","red")[factor(type)],
         col=color.palette,
         breaks = palette.breaks,
         scale = c("row"),
         dendrogram = "both")
dev.off()

################################

# 试一试cut tree，不同组细胞表达marker是不是有差别

#plot(hclust(RC.morula.corrected.top.variable))
complete.cluster = hclust(as.dist(1-abs(cor(log2(RC.morula.corrected.top.variable+1),method="spearman"))), 
                          method="complete")
plot(complete.cluster, hang = -1)

#reconcilePropertiesAndPrototype
rect.hclust(complete.cluster,3)


plot.exp = data.frame(exp = RC.morula.corrected.top.variable[which(rownames(RC.morula.corrected.top.variable) 
                                                                   == "BTBD10"),],
                      group = cutree(complete.cluster,3))

ggplot(plot.exp, aes(x = factor(group), y = log2(exp+1) )) +
  geom_boxplot()

