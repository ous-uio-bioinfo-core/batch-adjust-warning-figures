
# Figure 1


library(sva)
source("helperfunctions.r")
source("boxplot_function.r")
library(limma)
library(lsmeans)


set.seed(139) # selected based on illustrative values
ngenes = 1000 # Really not needed for removeBatchEffect, but a matrix is needed for ComBat.
index=1 # only this gene is plotted
usecolor=TRUE

sa = createsampleannotation(  list(c(25,5,0), c(0,20,50)), as.factors=TRUE) #  3 groups,  2 batches


matrix_true = matrix(rnorm(ngenes * nrow(sa), mean=10, sd=2), nrow=ngenes, ncol=nrow(sa))

# add a group effect
matrix_true[,sa$group==1] = matrix_true[,sa$group==1]+4
matrix_true[,sa$group==2] = matrix_true[,sa$group==2]+4
matrix_true[,sa$group==3] = matrix_true[,sa$group==3]+0


# calculating batchbox limits.
bbextra=1.2
batchboxheight=(max(matrix_true[index,]) - min(matrix_true[index,]) ) * bbextra
batchboxlow=min(matrix_true[index,]) - batchboxheight * (bbextra-1)/2
batchboxtop=max(matrix_true[index,]) + batchboxheight * (bbextra-1)/2
batchmeans = sapply( unique(sa$batch), FUN=function(x){  mean( matrix_true[index,sa$batch==x] ) } )
batchboxlowmeanoffsets=batchmeans - batchboxlow

# adding batch effects
matrix_batcheffect = matrix_true
matrix_batcheffect[,sa$batch==1] = matrix_batcheffect[,sa$batch==1]+2
matrix_batcheffect[,sa$batch==2] = matrix_batcheffect[,sa$batch==2]-2


# batch adjust with mean center.
matrix_batchadjusted = matrix_batcheffect
genemeans = rowMeans(matrix_batcheffect)
matrix_batchadjusted[,sa$batch ==1] =  genemeans +
  matrix_batchadjusted[,sa$batch ==1] - rowMeans( matrix_batchadjusted[,sa$batch ==1])
matrix_batchadjusted[,sa$batch ==2] =  genemeans + 
  matrix_batchadjusted[,sa$batch ==2] - rowMeans( matrix_batchadjusted[,sa$batch ==2])

# batch adjust with removeBatchEffect
mod = model.matrix(~as.factor(sa$group))
matrix_batchadjusted2 = removeBatchEffect(matrix_batcheffect[,], batch=as.factor(sa$batch), design=mod)


# plotting
figfilename = file.path( "../plots", paste("boxplots" , sep=""))
#figfile = paste( figfilename, ".pdf", sep=""); 
#pdf(file =figfile, width=8, height=4)

# eps
figfile = paste( figfilename, ".eps", sep="");
cairo_ps(file =figfile,  width=8, height=4)


# find the extreme ylims to be used by all plots.. 
tmpmatrix = rbind(matrix_true[index,], matrix_batcheffect[index,], matrix_batchadjusted[index,], matrix_batchadjusted2[index,])
ymin = min(t(t(tmpmatrix) - (matrix_true[index,]-batchboxlow)))
ymax = max(t(t(tmpmatrix) + (batchboxtop - matrix_true[index,])))
ylim = c(ymin, ymax)


layout(  matrix(c(1,2,3,4,5), 1, 5, byrow = TRUE) , widths=c(1,1,1,1,0.35) )
op=par(oma = c(5,4,0,0) + 0.1, mar = c(0,0,2,1) + 0.1, xpd=NA)
plot_one_gene(matrix_true[index,], group=sa$group,  batch=sa$batch,, 
              main=paste("(a) True values"), estimatemethod="CI", ylim=ylim,
              bbh=batchboxheight, bblo = batchboxlowmeanoffsets, usecolor=usecolor)
axis(side = 2, labels =TRUE)
plot_one_gene(matrix_batcheffect[index,], group=sa$group, batch=sa$batch,
              main=paste("(b) Batch affected values"), estimatemethod="CI", ylim=ylim,
              bbh=batchboxheight, bblo = batchboxlowmeanoffsets, usecolor=usecolor)
plot_one_gene(matrix_batchadjusted[index,], group=sa$group, batch=sa$batch,
              main=paste("(c) Mean centred values"), estimatemethod="CI", ylim=ylim,
              bbh=batchboxheight, bblo = batchboxlowmeanoffsets, usecolor=usecolor)
plot_one_gene(matrix_batchadjusted2[index,], group=sa$group, batch=sa$batch,
              main=paste("(d) Anova adjusted values"), estimatemethod="CI", ylim=ylim,
              bbh=batchboxheight, bblo = batchboxlowmeanoffsets, usecolor=usecolor)

estimatesboxesonly(matrix_batcheffect[index,], group=sa$group, batch=sa$batch,ylim=ylim, 
                   main="(e)", usecolor=usecolor)

dev.off()

print( paste("Figure created: ", normalizePath(figfile) ))

