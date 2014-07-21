
# Main script.
# Calls the individual scripts that creates figures 1,2,3a,3b.

source("http://bioconductor.org/biocLite.R")
biocLite( c("sva", "limma", "lsmeans", "bladderbatch"))


startwd = getwd()
source("overviewplot.r")

setwd("reanalysis/Leek2012")
source("figurescript_leek.r")
setwd(startwd)

setwd("reanalysis/Towfic2014")
source("figurescript_towfic.r")
setwd(startwd)

setwd("reanalysis/Johnson2007")
source("figurescript_johnson.r")
setwd(startwd)

sink("sessionInfo.txt")
sessionInfo()
sink()





