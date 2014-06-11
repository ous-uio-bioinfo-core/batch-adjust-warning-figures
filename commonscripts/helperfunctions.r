
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
adhocpvalueplot = function(realcombatp, reallimmap, randomp, main="P-values", xrange=1:25)
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


