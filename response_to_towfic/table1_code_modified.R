###########################################################################
#
#   Immuneering Corporation
#
#   SOFTWARE COPYRIGHT NOTICE AGREEMENT
#   This software and its documentation are copyright (2014) by the
#   Immuneering Corporation. All rights are reserved.
#
#   This software is supplied without any warranty or guaranteed support
#   whatsoever. Immuneering Corporation cannot be responsible for its use,
#   misuse, or functionality.
#
#
###########################################################################
#
#   The script has been modified from the original at
#      https://github.com/immuneering/biostat_reply
#   with alternative and corrected model specifications.
#
###########################################################################

library(utils)
library(lme4)
library(limma)
library(sva)
library(lmerTest) # for computing P values with lmer+anova
library(nlme) # lagt til for lme

# Download if not already exists
if (!exists('DataDir')) DataDir='.'
ensure.download <- function(url,destfile,mode='wb') {
  if (!file.exists(destfile)) {
    download.file(url=url,destfile=destfile,mode=mode)
  }  
}

getLabels <- function() {
  ensure.download(url='ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/GSE61901/matrix/GSE61901_series_matrix.txt.gz',
                       destfile = file.path(DataDir,'GSE61901_series_matrix.txt.gz'),mode = 'wb')
  data <- read.delim(file=file.path(DataDir,'GSE61901_series_matrix.txt.gz'),header=F,skip=32)
  SAMPLE <- as.character(unlist(data[2,-1]));
  CLASS <- as.character(unlist(data[26,-1]));
  SLOT <- as.character(unlist(data[12,-1]));
  SLOT <- gsub(pattern = 'array_address: ',replacement = '',x = SLOT)
  BATCH <- as.character(unlist(data[10,-1]));
  BATCH <- gsub(pattern = 'batch: ','B',BATCH)
  data.frame(SAMPLE=factor(SAMPLE),CLASS=factor(CLASS),SLOT=factor(SLOT),BATCH=factor(BATCH))
}

ensure.download(url = 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE61901&format=file&file=GSE61901%5Fnon%2Dnormalized%2Etxt%2Egz',
              destfile = file.path(DataDir,'GSE61901_non-normalized.txt.gz'),mode = 'wb')
data <- read.delim(file=file.path(DataDir,'GSE61901_non-normalized.txt.gz'),header=T,skip=5,row.names=1,as.is=T,check.names=F)
data <- data[,!(colnames(data) %in% c('Detection Pval'))]
data <- normalizeQuantiles(data)

# data=data[sample(1:nrow(data),100),] # Run to select probe subsample

labels <- getLabels()
labels <- labels[match(colnames(data),labels$SLOT),]
labels$chip <- factor(substr(as.character(labels$SLOT),1,nchar(as.character(labels$SLOT))-2))
combat_data <- ComBat(data,batch = labels$BATCH,mod = model.matrix(~labels$CLASS))
combat_data_averaged <- avearrays(combat_data,ID = labels$chip)
labels_averaged <- unique(labels[,c(2,4,5)]);
data_averaged <- avearrays(data,ID = labels$chip)

ranova_original <- function(X,lab,batch,sampid) {
  # Model used in letter by Towfic et al.
  # Random effects term incorrectly specified as Error(1/sampid)
  model <- aov(X~batch+lab+Error(1/sampid))
  unlist(summary(model))[['Error: Within.Pr(>F)2']]
}

ranova_equivalent <- function(X,lab,batch,sampid) {
  # Model equivalent to ranova_original: no error term
  model <- aov(X~batch+lab)
  unlist(summary(model))[['Pr(>F)2']]
}

ranova_corrected <- function(X,lab,batch,sampid) {
  # Correctly specified model with random effect term per sample
  # This becomes identical to twowayanova_blocking
  model <- aov(X~batch+lab+Error(sampid))
  unlist(summary(model))[['Error: sampid.Pr(>F)2']]
}

oneanova <- function(X,lab) {
  lab <- factor(as.character(lab))
  model <- glm(X~lab)
  coef(summary(model))[2,'Pr(>|t|)']
}

oneanova_alternative <- function(X,lab) {
  # Same as oneanove, but using aov instead of glm
  lab <- factor(as.character(lab))
  model <- aov(X~lab)
  unlist(summary(model))[['Pr(>F)1']]
}

twowayanova_blocking <- function(X,batch,lab) {
  batch <- factor(as.character(batch));
  lab <- factor(as.character(lab));
  model <- glm(X~batch+lab)
  coef(summary(model))[18,'Pr(>|t|)'] # NB: if number of batches changes, 18 must be changed
}

twowayanova_blocking_alternative <- function(X,batch,lab) {
  # Same as twowayanova_blocking, but using aov instead of glm
  # Extraction of P value from output does not depend on number of batches
  batch <- factor(as.character(batch));
  lab <- factor(as.character(lab));
  model <- aov(X~batch+lab)
  unlist(summary(model))[['Pr(>F)2']]
}

lmem_original <- function(X,lab,batch,sampid) {
  # Model used in letter by Towfic et al.
  # They made the following reference:
  #   Need to use ML estimates since REML estimates are not valid
  #   when fixed effects change
  #   http://stats.stackexchange.com/questions/41123/reml-vs-ml-stepaic
  # Direct model comparison of models with different fixed effects using REML
  # does not make sense. However, ML gives biased variance estimates. More
  # critically, the degrees of freedom does not take into account the mixed
  # random effects.
  model1 <- lmer(X~0+batch+(1|sampid),REML = F)
  model2 <- lmer(X~0+batch+lab+(1|sampid),REML=F);
  comp <- anova(model1,model2,test='Chisq');
  comp[['Pr(>Chisq)']][2]
}

lmem_corrected <- function(X,lab,batch,sampid) {
  # Corrected effect assessment using lmerTest
  # Should use REML for accurate variance estimation
  model <- lmer(X~0+batch+lab+(1|sampid),REML=T);
  anova(model)["lab","Pr(>F)"]
}

lmem_alternative <- function(X,lab,batch,sampid) {
  # Alternative model using nlme:lme
  # This fails on a few probes
  df = data.frame(X,batch,lab,sampid)
  tryCatch({
    model = lme(X ~ batch + lab, random= ~1|sampid, data=df)
    anova(model)['lab','p-value']
  },error=function(cond) {return(NA)})
}

limma_blocking <- function(inmat,lab) {
  class <- factor(make.names(as.character(lab$CLASS)))
  batch <- factor(make.names(as.character(lab$BATCH)))
  design <- model.matrix(~0+class+batch)
  fit <- lmFit(inmat,design)
  cm <- makeContrasts(DPvsQ=classGA.DP-classGA.Q,levels=design)
  fit2 <- contrasts.fit(fit, cm)
  limma.fit2 <- eBayes(fit2)
  limma.fit2$p.value
}

limma_noblocking <- function(inmat,lab) {
  class <- factor(make.names(as.character(lab$CLASS)))
  design <- model.matrix(~0+class)
  fit <- lmFit(inmat,design)
  cm <- makeContrasts(DPvsQ=classGA.DP-classGA.Q,levels=design)
  fit2 <- contrasts.fit(fit, cm)
  limma.fit2 <- eBayes(fit2)
  limma.fit2$p.value
}

limma_dupcorr <- function(inmat,lab) {
  batch <- factor(make.names(as.character(lab$BATCH)));
  class <- factor(make.names(as.character(lab$CLASS)))
  samearray <- lab$chip
  design <- model.matrix(~0+class+batch)
  corfit <- duplicateCorrelation(inmat,design,block=samearray)
  fit <- lmFit(inmat,design,block=samearray,correlation=corfit$consensus)
  cm <- makeContrasts(DPvsQ=classGA.DP-classGA.Q,levels=design)
  fit2 <- contrasts.fit(fit, cm)
  limma.fit2 <- eBayes(fit2)
  p.values <- limma.fit2$p.value
  attr(p.values,'consensus') <- corfit$consensus
  p.values
}

write.batch.sample.distribution <- function() {
  write.csv(as.matrix(table(unique_samples$labels.BATCH,unique_samples$labels.CLASS)),file='/workspace/copaxone/stage2/rna/ur_02112016_nygaard_data_balance/sample_table.csv',quote=F)
}

# Data with both strips kept as replicates
dp_vs_q <- list();
dp_vs_q$mat <- data[,labels$CLASS %in% c('GA DP','GA Q')]
dp_vs_q$labels <- labels[labels$CLASS %in% c('GA DP','GA Q'),]

# Data with slits averaged after ComBat normalisation
dp_vs_q_averaged <- list();
dp_vs_q_averaged$mat <- combat_data_averaged[,labels_averaged$CLASS %in% c('GA DP','GA Q')]
dp_vs_q_averaged$labels <- labels_averaged[labels_averaged$CLASS %in% c('GA DP','GA Q'),]

# Data with slits averaged
dp_vs_q_nocombat_averaged <- list();
dp_vs_q_nocombat_averaged$mat <- data_averaged[,labels_averaged$CLASS %in% c('GA DP','GA Q')]
dp_vs_q_nocombat_averaged$labels <- labels_averaged[labels_averaged$CLASS %in% c('GA DP','GA Q'),]


### Analyses that include strips are replicates
# Analyses from Towfic et al. letter
ranova_original_out <- apply(dp_vs_q$mat,1,ranova_original,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip)
lme4_original_out <- apply(dp_vs_q$mat,1,lmem_original,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip)
limma_out <- limma_dupcorr(dp_vs_q$mat,lab = dp_vs_q$labels)
# Equivalent model to ranova_original
ranova_equivalent_out <- apply(dp_vs_q$mat,1,ranova_equivalent,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip)
# Corrected models
ranova_corrected_out <- apply(dp_vs_q$mat,1,ranova_corrected,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip)
lme4_corrected_out <- apply(dp_vs_q$mat,1,lmem_corrected,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip)
# Mixed random effects model using nlme
lme4_alternative_out <- apply(dp_vs_q$mat,1,lmem_alternative,lab=dp_vs_q$labels$CLASS,batch=dp_vs_q$labels$BATCH,sampid=dp_vs_q$labels$chip)

### One-way models on ComBat batch-adjusted data
limma_noblocking_out <- limma_noblocking(inmat = dp_vs_q_averaged$mat,lab = dp_vs_q_averaged$labels)
oneanova_out <- apply(dp_vs_q_averaged$mat,1,oneanova,lab=dp_vs_q_averaged$labels$CLASS)

### Two-way models on non-adjusted data
limma_blocking_out <- limma_blocking(inmat = dp_vs_q_nocombat_averaged$mat,lab = dp_vs_q_nocombat_averaged$lab)
twowayanova_blocking_out <- apply(dp_vs_q_nocombat_averaged$mat,1,twowayanova_blocking,batch=dp_vs_q_nocombat_averaged$lab$BATCH,lab=dp_vs_q_nocombat_averaged$lab$CLASS)

### Two-way models on ComBat batch-adjusted data
limma_combat_blocking_out <- limma_blocking(inmat = dp_vs_q_averaged$mat,lab = dp_vs_q_nocombat_averaged$lab)
twowayanova_combat_blocking_out <- apply(dp_vs_q_averaged$mat,1,twowayanova_blocking,batch=dp_vs_q_nocombat_averaged$lab$BATCH,lab=dp_vs_q_nocombat_averaged$lab$CLASS)


# BH-correction and count number of predictions at FDR<limit
N.sign = function(P,method='BH',limit=0.05) if (sum(!is.na(P))) length(which(p.adjust(P,method=method) < limit)) else NA

result_summary <- data.frame(
  "Towfic et al.: technical replicates, CHIP blocking"=
    c(N.sign(limma_out),NA,N.sign(lme4_original_out),N.sign(ranova_original_out)),
  "Equivalent but simpler  model"=
    c(NA,NA,NA,N.sign(ranova_equivalent_out)),
  "Corrected models: Technical replicates, CHIP blocking"=
    c(NA,NA,N.sign(lme4_corrected_out),N.sign(ranova_corrected_out)),
  "Same model, different implementation"=
    c(NA,NA,N.sign(lme4_alternative_out),NA),
  "One-way: replicates averaged, ComBat normalized with treatment as covariate"=
    c(N.sign(limma_noblocking_out),N.sign(oneanova_out),NA,NA),
  "Two-way: replicates averaged, CHIP blocking"=
    c(N.sign(limma_blocking_out),N.sign(twowayanova_blocking_out),NA,NA),
  "Two-way: replicates averaged, ComBat normalized, CHIP blocking"=
    c(N.sign(limma_combat_blocking_out),N.sign(twowayanova_combat_blocking_out),NA,NA),
  "Result from Nygaard et al. article (GEO:GSEGSE40566)"=
    c(11,NA,NA,NA),
  row.names=c('limma','anova','mixed','aov+err'))

print(t(result_summary))
cat('Consensus correlation from duplicateCorrelation:',attr(limma_out,'consensus'))


# Data frame of all P values
P=data.frame(limma.dupcor=c(limma_out),ranova.orig=ranova_original_out,ranova.equiv=ranova_equivalent_out,lme4.orig=lme4_original_out,
             ranova=ranova_corrected_out,lme4=lme4_corrected_out,nlme=lme4_alternative_out,
             limma.oneway=c(limma_noblocking_out),anova.oneway=oneanova_out,
             limma.twoway=c(limma_blocking_out),anova.twoway=twowayanova_blocking_out,
             limma.combat=c(limma_combat_blocking_out),anova.combat=twowayanova_combat_blocking_out)

head(P) # P values for the first few probes

fmtprt <- function(tab) round(tab,7)
fmtprt(P[1:10,c(1,2,3,4)]) # ranova.orig = ranova.equiv
fmtprt(P[1:10,c(5,6,7,11)]) # ranova = anova.twoway ~ lme4 ~ nlme
fmtprt(P[1:10,c(8,9,10,12,13)])
