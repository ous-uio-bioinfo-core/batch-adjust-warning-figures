

library(sva)
library(bladderbatch)
data(bladderdata)
library(limma)

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

figfile = paste(getwd(), "/leekrandomdatapvalues.pdf", sep="")
pdf(file=figfile)
hist(pValuesComBat,  main="P-values, Random numbers", breaks=100, xlab="p-value")
dev.off()
print( paste("Figure created; ",figfile ))