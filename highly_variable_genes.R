library(genefilter)
library( statmod )




# 用之前要重新运行，以免被覆盖


##########

# Estimating technical noise

#RC.clean.clean.gene.DESeqN.morula 

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



# 
# # estimating the sample moments per gene
# meansMorula <- rowMeans( RC.clean.clean.gene.DESeqN.morula )
# varsMorula <- rowVars( RC.clean.clean.gene.DESeqN.morula )
# cv2Morula <- varsMorula / meansMorula^2
# 
# # exclude low expressing gene
# 
# minMeanForFit.Morula <- unname( quantile( meansMorula[ which( cv2Morula > .3 ) ], .95 ) )
# minMeanForFit.Morula
# 
# ## fit model
# 
# useForFit.Morula <- (meansMorula >= minMeanForFit.Morula)
# 
# fit.Morula <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansMorula[useForFit.Morula] ),
#                    cv2Morula[useForFit.Morula] )
# fit.Morula$coefficients
# 
# # actual noise coefficients 
# 
# 
# sfMorula = 
#   colData(dds.Gene)$sizeFactor[
#     match(colnames(RC.clean.clean.gene.DESeqN.morula),names(colData(dds.Gene)$sizeFactor))]
# 
# xi.Morula <- mean( 1 / sfMorula)
# 
# a0.morula <- unname( fit.Morula$coefficients["a0"] )
# a1.morula <- unname( fit.Morula$coefficients["a1tilde"] - xi.Morula )
# 
# 
# ## plot the fit
# colHeLa <- "#00207040"
# colAt <- "#70500040"
# colAtHi <- "#B0901040"
# 
# plot( NULL, xaxt="n", yaxt="n",
#       log="xy", xlim = c( 1e-1, 2e4 ), ylim = c( .005, 100 ),
#       xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)")
# 
# axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
#                              expression(10^4), expression(10^5) ) )
# 
# axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10","100" ), las=2 )
# 
# abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# 
# points( meansMorula, cv2Morula, pch=20, cex=.5, col=colHeLa )
# 
# 
# # Plot the fitted curve
# xg <- 10^seq( -2, 6, length.out=1000 )
# 
# lines( xg, (xi.Morula+a1.morula)/xg + a0.morula, col="#FF000080", lwd=3 )
# # Plot quantile lines around the fit
# 
# df <- ncol(RC.clean.clean.gene.DESeqN.morula) - 1
# 
# lines( xg, ( (xi.Morula+a1.morula)/xg + a0.morula  ) * qchisq( .975, df ) / df,
#        col="#FF000080", lwd=2, lty="dashed" )
# lines( xg, ( (xi.Morula+a1.morula)/xg + a0.morula  ) * qchisq( .025, df ) / df,
#        col="#FF000080", lwd=2, lty="dashed" )

#########
# perform the test for high variance genes
# 
# meansAt <- rowMeans( RC )
# varsAt <- rowVars( RC.clean.clean.ERCC.DESeqN.Morula )
# cv2At <- varsAt / meansAt^2
# 
# #match(colnames(RC.clean.clean.gene.DESeqN.morula),names(colData(dds.ERCC)$sizeFactor))
# 
# sfAt = colData(dds.ERCC)$sizeFactor[match(colnames(RC.clean.clean.gene.DESeqN.morula),names(colData(dds.ERCC)$sizeFactor))]
# 
# psia1theta <- mean( 1 / sfAt ) + a1.morula * mean( sfMorula / sfAt )
# 
# minBiolDisp <- .5^2
# 
# m <- ncol(RC.clean.clean.ERCC.DESeqN.Morula)
# 
# cv2th <- a0.morula + minBiolDisp + a0.morula * minBiolDisp
# 
# testDenom <- ( meansAt * psia1theta + meansAt^2 * cv2th ) / ( 1 + cv2th/m )
# 
# p <- 1 - pchisq( varsAt * (m-1) / testDenom, m-1 )
# 
# padj <- p.adjust( p, "BH" )
# sig <- padj < .1
# sig[is.na(sig)] <- FALSE
# table( sig )




##############################################
# perform actural test

meansAt <- rowMeans( RC.clean.clean.gene.DESeqN.morula )
varsAt <- rowVars( RC.clean.clean.gene.DESeqN.morula )
cv2At <- varsAt / meansAt^2
      
sfAt = 
   colData(dds.Gene)$sizeFactor[
     match(colnames(RC.clean.clean.gene.DESeqN.morula),names(colData(dds.Gene)$sizeFactor))]

psia1theta <- mean( 1 / sfAt ) + a1 * mean( sfHeLa / sfAt )

minBiolDisp <- .5^2

m <- ncol(RC.clean.clean.gene.DESeqN.morula)
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
########
      
sig.gene = names(sig)[sig]
write.table(sig.gene,file = "morula.high.variable.beforeCorrectEmbryoEffect.txt",quote = F,row.names = F)      

sig.big.gene = names(sig)[sig  & (rowMeans(RC.clean.clean.gene.DESeqN.morula) > 1)]
write.table(sig.big.gene,file = "morula.high.variable.bigthan5Readcount.beforeCorrectEmbryoEffect.txt",quote = F,row.names = F)      


