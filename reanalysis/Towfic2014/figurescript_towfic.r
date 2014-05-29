
# source("figuretowfic.r")
starttime = Sys.time()
debug = FALSE
downloaddata=TRUE
set.seed(100)


includelibs = c("Biobase", "GEOquery", "sva", "limma")
lapply(includelibs, require, character.only=T)
source("../../commonscripts/helperfunctions.r")
source("helperfunctions_towfic.r")

ret = loadtowfic(downloaddata)
sampleannotation = ret[["sampleannotation"]]
rawdata = ret[["data"]]

if(debug)
  rawdata = rawdata[1:1000,]

qnormdata = normalizeBetweenArrays(rawdata, method="quantile") 

# combat adjust
combatdata= as.matrix(ComBat(dat=qnormdata,
                             batch=sampleannotation$chip,
                             mod=model.matrix(~as.factor(sampleannotation$covariate)),
                             numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))

# Significance test DP vs N
group = factor(sampleannotation$covariate)
design = model.matrix(~0 + group)
fit = lmFit(combatdata, design)
cont.matrix = makeContrasts ( contrasts="groupDP-groupN", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
limma_p_alt = eBayes(fit2)$p.value[,1]
print(paste("ComBat adjusted real data, significant probes: ",
            sum(p.adjust(limma_p_alt, method="fdr")<0.05)))

#Limma blocked batch and significance test
group = factor(sampleannotation$covariate)
block = factor(sampleannotation$chip)
design = model.matrix(~0+group+block)
fit = lmFit(qnormdata, design)
cont.matrix = makeContrasts ( contrasts="groupDP-groupN", levels=design)  
fit2 = contrasts.fit(fit, cont.matrix)
limma_ret_woc = eBayes(fit2)
limma_p_woc = limma_ret_woc$p.value[,1]
print(paste("Limma adjusted real data, significant probes: ",  
            sum(p.adjust(limma_p_woc, method="fdr")<0.05)))

# radnom, ComBat adjusted
set.seed(100)
randdata = rawdata
randdata[,] = matrix(rnorm(length(randdata), mean=0, sd=1))
randcombatdata= as.matrix(ComBat(dat=randdata,
                                 batch=sampleannotation$chip,
                                 mod=model.matrix(~as.factor(sampleannotation$covariate)),
                                 numCovs=NULL, par.prior=TRUE, prior.plots=FALSE))
# Significance test DP vs N
group = factor(sampleannotation$covariate)
design = model.matrix(~0 + group)
fit = lmFit(randcombatdata, design)
cont.matrix = makeContrasts ( contrasts="groupDP-groupN", levels=design)
fit2 = contrasts.fit(fit, cont.matrix)
limma_p_rand_combat = eBayes(fit2)$p.value[,1]
print(paste("ComBat adjusted random data, significant probes: ",
            sum(p.adjust(limma_p_rand_combat, method="fdr")<0.05)))


# create pvalue plot
figfile = paste(getwd(), "/towficpvalues.pdf", sep="")
pdf(file=figfile)
adhocpvalueplot(limma_p_alt,limma_p_woc,limma_p_rand_combat)
dev.off()
print( paste("Figure created; ",figfile ))

print(paste( "Figure generated for Towfic et al data set. Time spent ", 
             as.integer(round(difftime(Sys.time(),starttime, units="mins"))  ), 
             " minutes", sep="") )

