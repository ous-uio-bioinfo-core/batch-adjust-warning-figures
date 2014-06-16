
# Creates a sampleannotation data.frame based on a list of batches with a vector of treatmentcounts.
# Used when creating artificial data.
createsampleannotation = function( treatment_in_batches)
{
	batches = vector()
	treatments = vector()
	for(i in 1:length(treatment_in_batches))
	{
		thisbatch = i
		for(s in 1:length(treatment_in_batches[[i]]))  
		{
			batches = c(batches, rep(thisbatch, treatment_in_batches[[i]][s]))
			treatments = c(treatments, rep( s, treatment_in_batches[[i]][s]))
		}
	}
	sampleannotation = data.frame(id=paste("sample", 1:length(treatments), sep=""),    treatment = treatments, batch=batches)
	return(sampleannotation)
}


# plots 3 histograms from p-value distibutions on top on each other.
# draws a line instead of bars for clarity.
# binsize is 0.01
oldadhocpvalueplot = function(realcombatp, reallimmap, randomp, main="P-values", xrange=1:25)
{
  thiscolors = c("red", "blue", "black")
  
  a = hist(realcombatp, breaks=100, plot=F)$counts[xrange]
  b = hist(reallimmap, breaks=100, plot=F)$counts[xrange]
  c = hist(randomp, breaks=100, plot=F)$counts[xrange]
  ylim=c(0, max(c(a,b,c)))
  # reproduced ComBat + limma
  plot((xrange)/100, a , ylim=ylim,
       main=main, xlab="p-value", ylab="frequency",
       type="l", lwd=2, col=thiscolors[1])
  
  # batch handled in limma
  lines((xrange)/100, b,
        col=thiscolors[2], lwd=2)
  
  # random ComBat + limma
  lines((xrange)/100, c,
        col=thiscolors[3], lwd=2)
  
  legend("topright", text.col=thiscolors,
         legend=c("Real data, ComBat adjusted",
                  "Real data, batch handled by Limma",
                  "Random data, ComBat adjusted")) 
}


# plots 3 histograms from p-value distibutions on top on each other.
# draws a line instead of bars for clarity.
# binsize is 0.01
adhocpvalueplot = function(realcombatp, reallimmap, randomp, main="P-values", xrange=1:25)
{

	a = hist(realcombatp, breaks=100, plot=F)$counts[xrange]
	b = hist(reallimmap, breaks=100, plot=F)$counts[xrange]
	c = hist(randomp, breaks=100, plot=F)$counts[xrange]
	ylim=c(0, max(c(a,b,c)))
	# reproduced ComBat + limma
	plot((xrange)/100, a , ylim=ylim,
			 main=main, xlab="p-value", ylab="frequency",
			 type="l", lwd=2, lty=1)
	
	# batch handled in limma
	lines((xrange)/100, b, lwd=2, lty=2)
	
	# random ComBat + limma
	lines((xrange)/100, c, lwd=2, lty=3)
	
	legend("topright", lwd=2, lty=1:3,
				 legend=c("Real data, ComBat adjusted",
				 				 "Real data, batch handled by Limma",
				 				 "Random data, ComBat adjusted")) 
}




# Run F test and store results in an object
Ftest = function (dat, mod, mod0,adjust=1)
{
	obj=new.env()
	class(obj)="Ftest"
	n <- dim(dat)[2]
	m <- dim(dat)[1]
	df1 <- dim(mod)[2]
	df0 <- dim(mod0)[2]
	p <- rep(0, m)
	Id <- diag(n)
	resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
	rss1 <- rowSums(resid * resid)
	rm(resid)
	resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*%  t(mod0))
	rss0 <- rowSums(resid0 * resid0)
	rm(resid0)
	fstats <- ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
	if (is.null(adjust)) {
		adjust=mean(rss1/(n - df1))/mean((rss0 - rss1)/(df1 - df0))
	}
	fadj=fstats*adjust
	p <- 1 - pf(fadj, df1 = (df1 - df0), df2 = (n - df1))
	obj$n1=n; obj$n0=m; obj$df1=df1-df0; obj$df0=n-df1;
	obj$F0=fstats; obj$F=fadj
	obj$MS1=(rss0-rss1)/(df1-df0)
	obj$MS0=rss1/(n-df1)
	obj$P=p
	obj$adjust=adjust
	return(obj)
}



