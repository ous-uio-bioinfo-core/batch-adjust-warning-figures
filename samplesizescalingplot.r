
### Figure of box-plots and confidence intervals

library(limma)
library(sva)

source("commonscripts/helperfunctions.r")

rseed=1111
set.seed(rseed)


pvalcounts = list()
batchaffectedgenefraction = 0.1
samplesizescaleing = c(1,10,100)
ngenes=20000

for(n in samplesizescaleing)
{
	sampleannotation = createsampleannotation(  list(c(n,n*5), c(n*5,n) )  ) 	
	
	x = table(sampleannotation[,2:3])	
	runname = paste(paste(x[,1], collapse="_"),  paste(x[,2], collapse="_"), sep="_x_")
	
	
	

	matrix_random = matrix(rnorm(ngenes * nrow(sampleannotation), mean=0, sd=1), nrow=ngenes, ncol=nrow(sampleannotation))
	for(batch in (unique(sampleannotation$batch)))
	{
    if((dim(matrix_random)[1]*batchaffectedgenefraction) >= 1)
    {
		  for(s in 1: (dim(matrix_random)[1]*batchaffectedgenefraction) )
		  {
		  	thisgenesbatcheffect = rnorm(1)
		  	a= sampleannotation$batch==batch
		  	matrix_random[s, a] = matrix_random[s, a] + thisgenesbatcheffect
		  }
    }
	}
	
	mod0 = model.matrix(~1,data=sampleannotation)
 	batch = factor(sampleannotation$batch)

  mod = model.matrix(~as.factor(group), data=sampleannotation)
	matrix_adjusted = ComBat(dat=matrix_random, batch=batch, 
 												mod=mod, numCovs=NULL, 
 												par.prior=TRUE, prior.plots=FALSE)
	pValues = f.pvalue(matrix_adjusted,mod,mod0)
	pvalcounts[["ComBat"]][[runname]]=hist(pValues, plot=FALSE, breaks=100)$counts
	
	# using removeBatchEffect
	#matrix_adjusted = removeBatchEffect(matrix_random, batch=batch, design=mod)
	#pValues = f.pvalue(matrix_adjusted,mod,mod0)
	#pvalcounts[["LIMMA"]][[runname]]=hist(pValues, plot=FALSE, breaks=100)$counts

}


for(i in 1:length(pvalcounts))
{
  figfilename = file.path( getwd(), "plots", paste("samplesizescaling",
                                                   names(pvalcounts)[i], 
                                                   sep="_"))
  figfile = paste( figfilename, ".pdf", sep=""); 
  pdf(file =figfile)
  
  
  
  plot(  (1:25)/100, 1:25, 
  		 main="(c) Effect of sample size scaling on P-values", 
  		 xlab="p-value", ylab="Frequency", type="n", 
  		 ylim=c(0, max(unlist( pvalcounts))))
  for(t in 1:length(pvalcounts[[i]]))
  {
  	lines((1:25)/100,pvalcounts[[i]][[t]][1:25],  lwd=2, lty=t)
  }
  
  legend("topright", legend=names(pvalcounts[[i]]), lty=1:length(pvalcounts[[i]]))
  
  
  dev.off()
  print( paste("Figure created: ",figfile ))
}
