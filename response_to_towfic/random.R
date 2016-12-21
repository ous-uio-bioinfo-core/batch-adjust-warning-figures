###
### Random data with mixed random effects
###

library(lme4)
library(lmerTest)
library(ggplot2)

## Specify matrix of N per group and batch, M replicates
# Towfic et al data: 2 groups, 17 batches
N.groups <- 2
N.group.batch <- matrix(c(4,0,0,2,1,1,5,0,4,0,1,1,4,1,2,1,2,1,1,0,1,2,1,2,2,0,2,0,2,0,1,0,1,0),N.groups)
M <- 2 # Number of replicates

N.batches <- ncol(N.group.batch)
N <- sum(N.group.batch) # number of samples
df <- data.frame(id=factor(rep(1:N,rep(M,N))))
df$group <- factor(rep(rep(1:2,N.batches),M*N.group.batch))
df$batch <- factor(rep(rep(1:N.batches,rep(2,N.batches)),M*N.group.batch))

sd.error <- .7 # Size of technical error term (relative to sample SD)
sd.batch <- .5 # Size of batch effect (relative to sample SD)

p.values <- data.frame()
for (i in 1:5000) {
  X.batch <- sd.batch*rnorm(N.batches)
  X.sample <- rnorm(N)
  df$x <- X.batch[df$batch] + X.sample[df$id] + sd.error*rnorm(N*M)
  # Likelihood ratio test (ML) as performed by Towfic et al
  # This will often produce small P values despite lack of a group difference!
  mod0 <- lme4::lmer(x ~ 0+batch+(1|id),df,REML=F)
  mod1 <- lme4::lmer(x ~ 0+batch+group+(1|id),df,REML=F)
  test.ml <- anova(mod1,mod0,test='chisq')
  p.ml <- test.ml[["Pr(>Chisq)"]][2]
  # F-test using approximated degrees of freedom
  mod <- lmerTest::lmer(x ~ batch+group+(1|id),df,REML=T)
  test.lmerTest <- anova(mod)
  p.lmerTest <- test.lmerTest[["Pr(>F)"]][2]
  # Add P values to table
  p.values <- rbind(p.values,data.frame(p.ml,p.lmerTest))
}
hist(p.values$p.ml,breaks=21)
hist(p.values$p.lmerTest,breaks=21)
#ggplot(p.values,aes(x=log10(p.lmerTest),y=log10(p.ml)))+geom_point()
