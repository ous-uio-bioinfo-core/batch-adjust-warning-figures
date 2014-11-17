
# Main script.
# Calls the individual scripts that creates figures 1,2,3a,3b.

# In order to reproduce the overview plot, version 3.20.8 of limma was needed.
# R version 3.1.1 was needed for installation of limma 3.20.8.
# The below two commands should download and install the nescesary packages if needed.
# source("http://bioconductor.org/biocLite.R")
# biocLite( c("sva", "limma", "lsmeans", "bladderbatch"))

# If R is not started with the appropriate folder as work directory:
#setwd("...path to batch-adjust-warning-figures folder...")

startwd = getwd()
source("overviewplot.r")

setwd("reanalysis/Leek2012")
source("figurescript_leek.r")
setwd(startwd)

source("samplesizescalingplot.r")

setwd("reanalysis/Towfic2014")
source("figurescript_towfic.r")
setwd(startwd)

setwd("reanalysis/Johnson2007")
source("figurescript_johnson.r")
setwd(startwd)

sink("sessionInfo.txt")
sessionInfo()
sink()





