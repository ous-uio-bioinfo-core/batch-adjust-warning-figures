
README for the batch-adjust-warning-figures git repository.
----------------------

This project consists of the r-code that produces the figures used in
["Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses"](http://biostatistics.oxfordjournals.org/content/early/2015/08/13/biostatistics.kxv027), 
Vegard Nygaard, Einar Andreas Rødland, and Eivind Hovig
Biostat first published online August 13, 2015
doi:10.1093/biostatistics/kxv027


The figures illustrate the side effects from batch effect adjustments using ComBat with study group as covariate or similar methods. 

### mainscript.r

Calls the individual figure scripts. The only script needed to be run in order to recreate the figures in the article. Takes about 6 minutes. 

### sessionInfo.txt

The output of sessionInfo() when mainscript was run. Be aware that the different versions of the r-packages could produce different results. In particular, the newest version of limma (3.20.8) was needed to produce the Figure 1.

### scripts/overviewplot.r

Figure 1. in the manuscript. 


### scripts/samplesizescalingplot.r

Figure 2c. in the manuscript. 

### scripts/helperfunctions.r boxplot_function.r

A few short utility functions.

### reanalysis/  

The individual data sets that are used to illustrate the adverse effect of ComBat.
Consist of a figure-script that loads the data, process the with or without ComBat, and runs a find differentially expressed genes analysis and creates a plot of the p-values. Figures 2a, 2b, 3a and 3c are created in their respective folder.


### plots/

The produced plots used in figure 1 and ficgure 2c.

### Additional resources

In our work, many more analyses than presented in the above figure-scripts were performed. These are more detailed and go beyond the scope of our article. But for the especially interested, some of them are made available as rmarkdown reports here:
https://github.com/ous-uio-bioinfo-core/batch-adjust-warning-reports.git
