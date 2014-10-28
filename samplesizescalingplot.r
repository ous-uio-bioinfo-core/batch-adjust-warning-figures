
### Figure of box-plots and confidence intervals

library(limma)
library(sva)

source("commonscripts/helperfunctions.r")

rseed=1111
set.seed(rseed)


pvalcounts = list()

for(n in c(1,10,100))
{
	sampleannotation = createsampleannotation(  list(c(n,n*5), c(n*5,n) )  ) 	
	
	x = table(sampleannotation[,2:3])	
	runname = paste(paste(x[,1], collapse="_"),  paste(x[,2], collapse="_"), sep="_x_")
	
	ngenes=20000
	batchaffectedgenefraction = 0.0001

	matrix_random = matrix(rnorm(ngenes * nrow(sampleannotation), mean=0, sd=1), nrow=ngenes, ncol=nrow(sampleannotation))
	for(batch in (unique(sampleannotation$batch)))
	{
		for(s in 1: (dim(matrix_random)[1]*batchaffectedgenefraction) )
		{
			thisgenesbatcheffect = rnorm(1)
			a= sampleannotation$batch==batch
			matrix_random[s, a] = matrix_random[s, a] + thisgenesbatcheffect
		}
	}
	
		mod0 = model.matrix(~1,data=sampleannotation)
 		batch = factor(sampleannotation$batch)
 		mod = model.matrix(~as.factor(treatment), data=sampleannotation)
# 	combat_edata = ComBat(dat=matrix_random, batch=batch, 
# 												mod=mod, numCovs=NULL, 
# 												par.prior=TRUE, prior.plots=FALSE)
	
	combat_edata = removeBatchEffect(matrix_random, batch=batch, design=mod)
	
	pValuesComBat = f.pvalue(combat_edata,mod,mod0)
	pvalcounts[[runname]]=hist(pValuesComBat, plot=FALSE, breaks=100)$counts
#qValuesComBat = p.adjust(pValuesComBat,method="BH")
}

figfilename = file.path( getwd(), "samplesizescaling")
figfile = paste( figfilename, ".pdf", sep=""); 
pdf(file =figfile, width=24, height=12)

#colorpal = c("blue", "red", "green", "black")
colorpal = c("black", "black","black","black")
plot(  (1:25)/100, 1:25, 
		 main="Effect of sample size scaling on P-values", 
		 xlab="p-value", ylab="Frequency", type="n", 
		 ylim=c(0, max(unlist( pvalcounts))))
for(i in 1:length(pvalcounts))
{
	lines((1:25)/100,pvalcounts[[i]][1:25], col=colorpal[i], lwd=2, lty=i)
}

legend("topright", legend=names(pvalcounts), text.col=colorpal, lty=1:length(pvalcounts))


dev.off()
print( paste("Figure created: ",figfile ))

