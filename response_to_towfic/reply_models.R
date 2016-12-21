###
### Analysis models
###

# Set to true to skip computationally demanding models
if (!exists('models.fast')) models.fast=FALSE

ranova <- function(X,lab,batch,sampid,type='original') {
  if (type=='original') {
    ### ORIGINAL CODE USED BY TOWFIC ET AL
    model <- aov(X~batch+lab+Error(1/sampid)) # no effect of Error-term
    unlist(summary(model))[['Error: Within.Pr(>F)2']]
  } else if (type=='equivalent') {
    ### EQUIVALENT MODEL
    model <- aov(X~batch+lab) # actual model when using misspecified error term
    unlist(summary(model))[['Pr(>F)2']]
  } else if (type=='correct') {
    ### CORRECTED MODEL
    model <- aov(X~batch+lab+Error(sampid)) # becomes identical to twowayanova_blocking
    unlist(summary(model))[['Error: sampid.Pr(>F)2']]
  } else NA
}

oneanova <- function(X,lab,type='aov') {
  lab <- factor(as.character(lab))
  if (type=='glm') {
    model <- glm(X~lab)
    coef(summary(model))[2,'Pr(>|t|)']
  } else if (type=='aov') {
    model <- aov(X~lab)
    unlist(summary(model))[['Pr(>F)1']]
  } else NA
}

twowayanova_blocking <- function(X,batch,lab,type='aov') {
  batch <- factor(as.character(batch));
  lab <- factor(as.character(lab));
  if (type=='glm') { # NB: index 18 depends on number of batches in analysis!
    model <- glm(X~batch+lab)
    coef(summary(model))[18,'Pr(>|t|)']
  } else if (type=='aov') { # Same model using aov
    model <- aov(X~batch+lab)
    unlist(summary(model))[['Pr(>F)2']]
  } else NA
}

lmem <- function(X,lab,batch,sampid,type='original') {
  if (models.fast) {
    NA
  } else if (type=='original') {
    # Need to use ML estimates since REML estimates are not valid
    # when fixed effects change
    # http://stats.stackexchange.com/questions/41123/reml-vs-ml-stepaic
    model1 <- lme4::lmer(X~0+batch+(1|sampid),REML = F)
    model2 <- lme4::lmer(X~0+batch+lab+(1|sampid),REML=F);
    comp <- anova(model1,model2,test='Chisq');
    comp[['Pr(>Chisq)']][2]
  } else if (type=='correct') {
    # Corrected effect assessment using lmerTest
    # Should used REML for accurate variance estimation
    model <- lmerTest::lmer(X~0+batch+lab+(1|sampid),REML=T);
    anova(model)["lab","Pr(>F)"]
  } else if (type=='nlme') {
    # Alternative model using nlme:lme
    df = data.frame(X,batch,lab,sampid)
    tryCatch({
      model = lme(X ~ batch + lab, random= ~1|sampid, data=df)
      anova(model)['lab','p-value']
    },error=function(cond) {return(NA)})
  } else NA
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
  cat("Consensus correlation: ",corfit$consensus)
  fit <- lmFit(inmat,design,block=samearray,correlation=corfit$consensus)
  cm <- makeContrasts(DPvsQ=classGA.DP-classGA.Q,levels=design)
  fit2 <- contrasts.fit(fit, cm)
  limma.fit2 <- eBayes(fit2)
  P.vals=limma.fit2$p.value
  attr(P.vals,'consensus') <- corfit$consensus
  return(P.vals)
}

