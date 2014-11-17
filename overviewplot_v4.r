


library(sva)
source("commonscripts/helperfunctions.r")
source("commonscripts/boxplot_function_v4.r")
library(limma)
library(lsmeans)






#set.seed(100)
ngenes = 1000
index=1
#sa = createsampleannotation(  list(c(10,5,0,0), c(0,5,5,0), c(0,0,5,10)), as.factors=TRUE)


sa = createsampleannotation(  list(c(10,3,0,0), c(0,10,30,0), c(0,0,3,20)), as.factors=TRUE)


#sa = createsampleannotation(  list(c(20,2,0), c(0,20,50)), as.factors=TRUE) # god 3 grupper 2 batcher


#sa = createsampleannotation(  list(c(15,0), c(2,2), c(0,15)), as.factors=TRUE)
matrix_true = matrix(rnorm(ngenes * nrow(sa), mean=5, sd=2), nrow=ngenes, ncol=nrow(sa))

# add a group effect for group 3
 matrix_true[,sa$group==1] = matrix_true[,sa$group==1]+4
 matrix_true[,sa$group==2] = matrix_true[,sa$group==2]+4
 matrix_true[,sa$group==3] = matrix_true[,sa$group==3]+0
 matrix_true[,sa$group==4] = matrix_true[,sa$group==4]+0

bbextra=1.2
batchboxheight=(max(matrix_true[index,]) - min(matrix_true[index,]) ) * bbextra
batchboxlow=min(matrix_true[index,]) - batchboxheight * (bbextra-1)/2
batchboxtop=max(matrix_true[index,]) + batchboxheight * (bbextra-1)/2
batchmeans = sapply( unique(sa$batch), FUN=function(x){  mean( matrix_true[index,sa$batch==x] ) } )
batchboxlowmeanoffsets=batchmeans - batchboxlow


matrix_batcheffect = matrix_true
matrix_batcheffect[,sa$batch==1] = matrix_batcheffect[,sa$batch==1]+2
matrix_batcheffect[,sa$batch==2] = matrix_batcheffect[,sa$batch==2]-2


#Batch adjust with mean center.

matrix_batchadjusted = matrix_batcheffect
genemeans = rowMeans(matrix_batcheffect)
matrix_batchadjusted[,sa$batch ==1] =  genemeans +
  matrix_batchadjusted[,sa$batch ==1] - rowMeans( matrix_batchadjusted[,sa$batch ==1])
matrix_batchadjusted[,sa$batch ==2] =  genemeans + 
  matrix_batchadjusted[,sa$batch ==2] - rowMeans( matrix_batchadjusted[,sa$batch ==2])


mod = model.matrix(~as.factor(sa$group))
matrix_batchadjusted2 = removeBatchEffect(matrix_batcheffect[,], batch=as.factor(sa$batch), design=mod)

source("commonscripts/boxplot_function_v4.r")
figfilename = file.path( getwd(), "plots", paste("boxplots_v4" , sep=""))
figfile = paste( figfilename, ".pdf", sep=""); 
#pdf(file =figfile, width=6, height=6)

# find the ylim.
tmpmatrix = rbind(matrix_true[index,], matrix_batcheffect[index,], matrix_batchadjusted[index,], matrix_batchadjusted2[index,])
ymin = min(t(t(tmpmatrix) - (matrix_true[index,]-batchboxlow)))
ymax = max(t(t(tmpmatrix) + (batchboxtop - matrix_true[index,])))
ylim = c(ymin, ymax)




#op=par(mfrow=c(3, 2),   xaxt="n", mar=c(1,3,2,1)+0.1)
layout(  matrix(c(1,2,0,3,4,5), 2, 3, byrow = TRUE) , widths=c(1,1,0.35) )
layout.show(5)
op=par(oma = c(5,4,0,0) + 0.1, mar = c(0,0,2,1) + 0.1)
plot_one_gene(matrix_true[index,], group=sa$group,  batch=sa$batch,, 
              main=paste("(a) True values"), estimatemethod="CI", ylim=ylim,
              bbh=batchboxheight, bblo = batchboxlowmeanoffsets)
axis(side = 2, labels =TRUE)
plot_one_gene(matrix_batcheffect[index,], group=sa$group, batch=sa$batch,
              main=paste("(b) Batch affected values"), estimatemethod="CI", ylim=ylim,
              bbh=batchboxheight, bblo = batchboxlowmeanoffsets)
plot_one_gene(matrix_batchadjusted[index,], group=sa$group, batch=sa$batch,
              main=paste("(c) Mean adjusted values"), estimatemethod="CI", ylim=ylim,
              bbh=batchboxheight, bblo = batchboxlowmeanoffsets)
axis(side = 2, labels =TRUE)
plot_one_gene(matrix_batchadjusted2[index,], group=sa$group, batch=sa$batch,
              main=paste("(d) Anova adjusted values"), estimatemethod="CI", ylim=ylim,
              bbh=batchboxheight, bblo = batchboxlowmeanoffsets)

#plot_one_gene(matrix_batchadjusted2[index,], group=sa$group, batch=sa$batch, ylim=ylim,
#              main=paste("(e) Batch affected values, batch included in ANOVA"),estimatemethod="lsmeans",
#              bbh=batchboxheight, bblo = batchboxlowmeanoffsets)
estimatesboxesonly(matrix_batcheffect[index,], group=sa$group, batch=sa$batch,ylim=ylim, 
                   main="(e) 2-Way-Anova")



## show the regions that have been allocated to each plot




#dev.off()
