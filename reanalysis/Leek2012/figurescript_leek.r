

# Script to create a Figure 2a and 2b.

# This is a modified version of the ComBat use example from the sva tutorial:
# http://www.bioconductor.org/packages/2.13/bioc/vignettes/sva/inst/doc/sva.pdf
# The modifications consist of swapping the data with random numbers a plot of the p-value distibution and qq-plot of the F-statisitcs

library(sva)
library(bladderbatch)
data(bladderdata)
library(limma)
source("../../scripts/helperfunctions.r")

pheno = pData(bladderEset)
edata = exprs(bladderEset)

# Substitute the data with random numbers from a normal distribution.
set.seed(100)
edata[,] = rnorm(length(edata), mean=0, sd=1)


mod0 = model.matrix(~1,data=pheno)

batch = pheno$batch

mod = model.matrix(~as.factor(cancer), data=pheno)

combat_edata = ComBat(dat=edata, batch=batch, 
                      mod=mod, numCovs=NULL, 
                      par.prior=TRUE, prior.plots=FALSE)

pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

#figfile = paste(getwd(), "/leekrandomdatapvalues.pdf", sep="")
#pdf(file=figfile)

# eps
figfile = paste(getwd(), "/leekrandomdatapvalues.eps", sep="")
cairo_ps(file =figfile)


hist(pValuesComBat,  main="(a) P-values, Random numbers", breaks=100, xlab="p-value")
dev.off()
print( paste("Figure created; ",figfile ))



#figfile = paste(getwd(), "/leekqqplot.png", sep="")
#png(file=figfile)

# eps
figfile = paste(getwd(), "/leekqqplot.eps", sep="")
cairo_ps(file =figfile)

test = Ftest(combat_edata,mod,mod0)
Fquant=qf(ppoints(nrow(combat_edata)),test$df1,test$df0)
qqplot(Fquant,test$F, ylim=c(0,max(test$F)), xlim=c(0,max(test$F)),
       xlab='Quantiles of F distribution',
       ylab='F statistics from analyses',
       main="(b) QQ-plot of F-statistics")
abline(0,1,lwd=2, lty=2)
dev.off()
print( paste("Figure created; ", normalizePath(figfile) ))
