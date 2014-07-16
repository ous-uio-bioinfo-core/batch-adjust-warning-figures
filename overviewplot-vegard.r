

# Script to create Figure 1

library(sva)
library(lsmeans)

adhocboxplot = function(y, grouplabels, batchlabels, addlegend=FALSE, ylim=NULL, xlab="x", figureletter=NA)
{
  #x = batchlabels + grouplabels/8
  colpalette = c("blue", "red", "darkgreen", "magenta")
  pchselection = c(1, 2, 3)
  x = batchlabels
  lwd=2
  
  # add extra x for better visual separation
  for(b in unique(batchlabels))
  {
    
    separation = 0.15
    xoffset=0    
    for(g in unique(grouplabels[batchlabels==b]))
    {
      x[batchlabels==b & grouplabels==g] = x[batchlabels==b & grouplabels==g] + xoffset
      xoffset = xoffset+separation
    }
  }
  
  xlim = c(0, length(unique(batchlabels))+2)
  plot(x, y, col=colpalette[grouplabels], pch=pchselection[batchlabels], 
       xlim=xlim, lwd=lwd, cex=2, ylim=ylim, xlab=xlab, ylab="", axes=FALSE, cex.lab =2)
  
  #add boxplot
  groupnames = unique(grouplabels)
  groupnames = groupnames[order(groupnames)]
  xoffset=length(unique(batchlabels))+1
  boxwidth=0.35
  boxseparation = 0.25
  for(g in groupnames)
  { 

    m = mean(y[grouplabels==g])
    ci = t.test(y[grouplabels==g])[["conf.int"]]
    xleft = xoffset
    ybottom = ci[1]
    xright = xleft + boxwidth
    ytop = ci[2]
    rect(xleft, ybottom, xright, ytop, border=colpalette[g], lwd=lwd)
    lines(c(xoffset,xright), c(m,m), col=colpalette[g], lwd=lwd)
    xoffset = xoffset+boxseparation
  }
  legendcex=2
  if(addlegend)
  {
    
    legend("topright", legend=paste("Group ", groupnames, sep=""), 
           text.col=colpalette[groupnames], bty="n", cex=legendcex)
    legend("bottomright", legend=paste("Batch ", unique(batchlabels), sep=""),
           pch=pchselection[unique(batchlabels)], bty="n", cex=legendcex)
  }
  if(!is.na(figureletter))
  {
  	figureletter=paste("(", figureletter, ")", sep="")
  	text(x=0, y=ylim[2], labels=figureletter, cex=legendcex)
  }
}



adhocboxplot2 = function(y, grouplabels, batchlabels, addlegend=FALSE, ylim=NULL, xlab="x", figureletter=NA)
{
	#x = batchlabels + grouplabels/8
	colpalette = c("blue", "red", "darkgreen", "magenta")
	pchselection = c(1, 2, 3)
	x = batchlabels
	lwd=2
	
	
	xlim = c(0, length(unique(batchlabels))+2)
	plot(1, 1, col=colpalette[grouplabels], pch=pchselection[batchlabels], 
			 xlim=xlim, lwd=lwd, cex=2, ylim=ylim, xlab=xlab, ylab="", axes=FALSE, cex.lab =2, type="n")
	
	expression =  matrix_conditionbatch[1,]
	group= factor(grouplabels)
	batch= factor(batchlabels)
	fit = lm(expression ~ group + batch )
	anovaest = summary(lsmeans(fit,  ~group))
	
	
	groupnames = unique(grouplabels)
	groupnames = groupnames[order(groupnames)]
	xoffset=2
	boxwidth=0.35
	boxseparation = 0.25
	for(g in groupnames)
	{ 
		m =anovaest[g , "lsmean"]		
		xleft = xoffset
		ybottom = anovaest[g , "lower.CL"]		
		xright = xleft + boxwidth
		ytop = anovaest[g , "upper.CL"]	
		rect(xleft, ybottom, xright, ytop, border=colpalette[g], lwd=lwd)
		lines(c(xoffset,xright), c(m,m), col=colpalette[g], lwd=lwd)
		xoffset = xoffset+boxseparation
	}
	legendcex=2
	if(addlegend)
	{
		
		legend("topright", legend=paste("Group ", groupnames, sep=""), 
					 text.col=colpalette[groupnames], bty="n", cex=legendcex)
		legend("bottomright", legend=paste("Batch ", unique(batchlabels), sep=""),
					 pch=pchselection[unique(batchlabels)], bty="n", cex=legendcex)
	}
	if(!is.na(figureletter))
	{
		figureletter=paste("(", figureletter, ")", sep="")
		text(x=0, y=ylim[2], labels=figureletter, cex=legendcex)
	}
}


source("commonscripts/helperfunctions.r")

# 3 batches.
# 4 groups
#sampleannotation = createsampleannotation(  list(c(100,0,0,20), c(0,100,0,20), c(0,0,100,20))) 
sampleannotation = createsampleannotation(  list(c(50,0,0,20), c(0,50,0,20), c(0,0,50,20))) 

table(sampleannotation[,2:3])

ngenes=1000
#matrix_random = createrandomdata(ngenes, sampleannotation, mean=0, sd=1)
set.seed(6)
matrix_random = matrix(rnorm(ngenes * nrow(sampleannotation), mean=0, sd=1), nrow=ngenes, ncol=nrow(sampleannotation))

matrix_condition = matrix_random
matrix_condition[, sampleannotation$treatment %in% c(3,4)] = matrix_condition[, sampleannotation$treatment %in% c(3,4)] +2
matrix_conditionbatch=matrix_condition
matrix_conditionbatch[, sampleannotation$batch ==2] = matrix_condition[, sampleannotation$batch ==2] -1
matrix_conditionbatch[, sampleannotation$batch ==3] = matrix_condition[, sampleannotation$batch ==3] +1
matrix_meancenter = matrix_conditionbatch
matrix_meancenter[,sampleannotation$batch ==1] =  matrix_meancenter[,sampleannotation$batch ==1] - rowMeans( matrix_meancenter[,sampleannotation$batch ==1])
matrix_meancenter[,sampleannotation$batch ==2] =  matrix_meancenter[,sampleannotation$batch ==2] - rowMeans( matrix_meancenter[,sampleannotation$batch ==2])
matrix_meancenter[,sampleannotation$batch ==3] =  matrix_meancenter[,sampleannotation$batch ==3] - rowMeans( matrix_meancenter[,sampleannotation$batch ==3])
mod = model.matrix(~as.factor(sampleannotation$treatment))
matrix_combat = ComBat(dat=matrix_conditionbatch, batch=sampleannotation$batch,
                       mod=mod, par.prior=TRUE, prior.plots=FALSE)




figfile = paste( getwd(), "/", "boxplots",  ".pdf", sep="")
pdf(file =figfile, width=32, height=16)
par(mfrow=c(1, 5))

alldata = c(matrix_condition[index,],
            matrix_conditionbatch[index,],
            matrix_meancenter[index,],
            matrix_combat[index,])
ylim = c(min(alldata), max(alldata))
adhocboxplot(matrix_condition[index,], sampleannotation$treatment, sampleannotation$batch, addlegend=FALSE, ylim=ylim, xlab="Without batch effects", figureletter="a")
adhocboxplot(matrix_conditionbatch[index,], sampleannotation$treatment, sampleannotation$batch, addlegend=FALSE, ylim=ylim, xlab="With batch effects added", figureletter="b")
adhocboxplot(matrix_meancenter[index,], sampleannotation$treatment, sampleannotation$batch, addlegend=FALSE, ylim=ylim, xlab="Zero-centered per batch", figureletter="c")
adhocboxplot(matrix_combat[index,], sampleannotation$treatment, sampleannotation$batch, addlegend=FALSE, ylim=ylim, xlab="ANOVA centered values (for tiden combat)", figureletter="d")
adhocboxplot2(matrix_conditionbatch[index,], sampleannotation$treatment, sampleannotation$batch, addlegend=TRUE, ylim=ylim, xlab="ANOVA estimates", figureletter="e")
dev.off()

print( paste("Figure created; ",figfile ))



