
## May wish to set working directory and/or datadir
# setwd(...) # Working directory containing scripts
# DataDir=... # Where data is store permanently

library(utils)
library(sva)
library(limma)
library(lme4) # original lmer
library(lmerTest) # used to get P values using lmer+anova
library(nlme) # used for lme

# Include scripts
source('reply_models.R')
source('reply_subroutines.R')
source('reply_data.R')


### Data selection
#NAME='All'; models.fast=F; # All methods included
NAME='Some'; models.fast=T; # Only some of the methods included
#plot(data.all.avg,data.all.sd,cex=0.1)
data <- select_data(1000) # sample=N to subsample
#data <- select_data(sd.min=.5) # select subset av expression avg or sd
#data <- select_data(subset=which(data.all.sd<.5),subset.tags="Slt0.5")
cat("Data size:",dim(data))

# Get labels for selected data
labels <- labels.all
labels <- labels[match(colnames(data),labels$SLOT),]
labels$chip <- factor(substr(as.character(labels$SLOT),1,nchar(as.character(labels$SLOT))-2))

labels$include <- labels$CLASS %in% c('GA DP','GA Q')
if (F) { # only include batches that contain samples from both groups
  commonbatches <- intersect(labels$BATCH[labels$CLASS=='GA DP'],labels$BATCH[labels$CLASS=='GA Q'])
  labels$include <- labels$include & labels$BATCH %in% commonbatches
  attr(data,'log')$tags <- c(attr(data,'log')$tags,'commonbatches') 
}

# Randomise data
#data <- permute_dp_q_samples(data,labels,rows=T,seed=1004)
Log.print(data)


### Continue normalisation and averaging
combat_data <- ComBat(data,batch = labels$BATCH,mod = model.matrix(~labels$CLASS))
combat_data_averaged <- avearrays(combat_data,ID = labels$chip)
labels_averaged <- unique(labels[,c(2,4,5,6)]);
data_averaged <- avearrays(data,ID = labels$chip)


## Final data sets for analysis
dp_vs_q <- list();
dp_vs_q$mat <- data[,labels$include]
dp_vs_q$labels <- labels[labels$include,]

dp_vs_q_averaged <- list();
dp_vs_q_averaged$mat <- combat_data_averaged[,labels_averaged$include]
dp_vs_q_averaged$labels <- labels_averaged[labels_averaged$include,]

dp_vs_q_nocombat_averaged <- list();
dp_vs_q_nocombat_averaged$mat <- data_averaged[,labels_averaged$include]
dp_vs_q_nocombat_averaged$labels <- labels_averaged[labels_averaged$include,]

# Original models
ranova_orig_out <- apply(dp_vs_q$mat,1,ranova,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip,type='original')
ranova_equiv_out <- apply(dp_vs_q$mat,1,ranova,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip,type='equivalent')
lme4_orig_out <- apply(dp_vs_q$mat,1,lmem,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip,type='original')

# Adjusted/corrected models
ranova_out <- apply(dp_vs_q$mat,1,ranova,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip,type='correct') # type=2 is correct model
lme4_out <- apply(dp_vs_q$mat,1,lmem,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip,type='correct') # type=1 is correct model
nlme_out <- apply(dp_vs_q$mat,1,lmem,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip,type='nlme') # type=2 is nlme model
limma_out <- limma_dupcorr(dp_vs_q$mat,lab = dp_vs_q$labels)

# 1-way ANOVA
limma_noblocking_out <- limma_noblocking(inmat = dp_vs_q_averaged$mat,lab = dp_vs_q_averaged$labels)
oneanova_out <- apply(dp_vs_q_averaged$mat,1,oneanova,lab=dp_vs_q_averaged$labels$CLASS)

# 2-way ANOVA
limma_blocking_out <- limma_blocking(inmat = dp_vs_q_nocombat_averaged$mat,lab = dp_vs_q_nocombat_averaged$labels)
twowayanova_blocking_out <- apply(dp_vs_q_nocombat_averaged$mat,1,twowayanova_blocking,batch=dp_vs_q_nocombat_averaged$labels$BATCH,lab=dp_vs_q_nocombat_averaged$labels$CLASS)

# 2-way ANOVA after ComBat
limma_combat_blocking_out <- limma_blocking(inmat = dp_vs_q_averaged$mat,lab = dp_vs_q_averaged$labels)
twowayanova_combat_blocking_out <- apply(dp_vs_q_averaged$mat,1,twowayanova_blocking,batch=dp_vs_q_averaged$labels$BATCH,lab=dp_vs_q_averaged$labels$CLASS)


table(factor(dp_vs_q_averaged$labels$CLASS))
table(factor(dp_vs_q_averaged$labels$BATCH),factor(dp_vs_q_averaged$labels$CLASS))

result_summary <- data.frame(
  "Original, faulty, models with replicates and blocking"=
    c(N.sign(limma_out),N.sign(ranova_orig_out),N.sign(lme4_orig_out),NA),
  "Actual model, equivalent to original model"=
    c(NA,N.sign(ranova_equiv_out),NA,NA),
  "Corrected model with replicates and CHIP blocking"=
    c(NA,N.sign(ranova_out),N.sign(lme4_out),N.sign(nlme_out)),
  "One-way: averaged replicates, ComBat with treatments as covariates"=
    c(N.sign(limma_noblocking_out),N.sign(oneanova_out),NA,NA),
  "Two-way: averaged replicates, CHIP blocking"=
    c(N.sign(limma_blocking_out),N.sign(twowayanova_blocking_out),NA,NA),
  "Two-way: averaged replicates, ComBat, CHIP blocking"=
    c(N.sign(limma_combat_blocking_out),N.sign(twowayanova_combat_blocking_out),NA,NA),
  row.names = c('limma','anova','lme4','nlme')
)
print(t(result_summary))
cat('Consensus correlation:',attr(limma_out,'consensus'))

### Similarity between methods
P=data.frame(limma.dupcor=c(limma_out),ranova.orig=ranova_orig_out,ranova.actual=ranova_equiv_out,lme4.orig=lme4_orig_out,
             ranova=ranova_out,lme4=lme4_out,nlme=nlme_out,
             limma.oneway=c(limma_noblocking_out),anova.oneway=oneanova_out,
             limma.twoway=c(limma_blocking_out),anova.twoway=twowayanova_blocking_out,
             limma.combat=c(limma_combat_blocking_out),anova.combat=twowayanova_combat_blocking_out)
P.nonmissing <- apply(!is.na(P),2,sum)>0
P.red <- P[,P.nonmissing]
P.dist=dist2(-log(1e-10+P.red))
-log(P.dist)/log(10)*10

plot.sample=sample(nrow(P.red),min(nrow(P.red),1000))
plot.sample=which(apply(data,1,sd)>.5)
pairs((P.red[plot.sample,]),diag.panel=panel.hist,cex=.1)

if (FALSE) { # coloured subgroups
  plot.grp=apply(data,1,mean)>7
  pairs((P.red[plot.sample,]),diag.panel=panel.hist,bg=c('red','green')[factor(plot.grp[plot.sample])],cex=.1,pch=21)
}

### Save results
#store.results(get.tag.file(data,name=NAME))
