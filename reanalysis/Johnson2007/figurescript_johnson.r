
# source("figurescriptjohnson.r")
starttime = Sys.time()

library("sva")
library("limma")
source("../../commonscripts/helperfunctions.r")
source("helperfunctions_johnson.r")

tmp = loadjohnsondata()
sampleannotation = tmp[["sampleannotation"]]
normdata = tmp[["normdata"]]

combatdata = as.matrix( ComBat(
  dat=normdata, 
  batch=sampleannotation$Batch, 
  mod=model.matrix(~as.factor(sampleannotation$"Type")), 
  numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))

# DE pvalues ComBat adjusted real data
Type = as.factor(sampleannotation$Type)
design = model.matrix(~0 + Type)
fit = lmFit(combatdata, design)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
combatp = eBayes(fit2)$p.value[,1]
rm(Type, design, fit, cont.matrix, fit2)
print(paste("ComBat adjusted real data, significant probes: ",  sum(p.adjust(combatp, "fdr")<0.05)))

# DE pvalues not ComBat adjusted, 
# batch blocked in limma
Type = as.factor(sampleannotation$Type)
Block = as.factor(sampleannotation$Batch)
design = model.matrix(~0+Type+Block)
fit = lmFit(normdata, design)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
limmap = eBayes(fit2)$p.value[,1]
rm(Type, Block, design, fit, cont.matrix, fit2)
print(paste("Limma batch adjusted real data, significant probes: ",  sum(p.adjust(limmap, "fdr")<0.05)))

#Random data
set.seed(100)
randdata = normdata
randdata[,] =rnorm(length(randdata), mean=0, sd=1)
randcombatdata = as.matrix(ComBat(
  dat=randdata, 
  batch=sampleannotation$Batch, 
  mod=model.matrix(~as.factor(sampleannotation$Type)), 
  numCovs=NULL, par.prior=TRUE,  prior.plots=FALSE))

Type = as.factor(sampleannotation$Type)
design = model.matrix(~0 + Type)
cont.matrix = makeContrasts ( contrasts="TypeR-TypeC", levels=design) 
fit = lmFit(randcombatdata, design)
fit2 = contrasts.fit(fit, cont.matrix)
randp = eBayes(fit2)$p.value[,1]
print(paste("ComBat adjusted random data, significant probes: ",  sum(p.adjust(randp, "fdr")<0.05)))

# create pvalue plot
figfile = paste(getwd(), "/dataset2pvalues.pdf", sep="")
pdf(figfile)
adhocpvalueplot( combatp, limmap, randp,  main="(b) P-values")
dev.off()
print( paste("Figure created; ",figfile ))

print(paste( "Figure generated for pvalues from Johnson data set 2. Time spent ", 
             as.integer(round(difftime(Sys.time(),starttime, units="mins"))  ), 
             " minutes", sep="") )





      