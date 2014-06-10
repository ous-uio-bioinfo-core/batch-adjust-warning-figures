
README for the batch adjust warning figures git repository.
----------------------

This project consists of the r-code that produces the figures used in "Methods that remove batch effects while retaining group
differences may lead to exaggerated confidence in downstream analyses", V. Nygaard, E. A.
Rødland, E. Hovig, manuscript in preparation.

The figures illustrate the side effects from batch effect adjustments using ComBat with study group as covariate or similar methods.

### figurescript.r

Calls the individual figure scripts. The only script needed to be run in order to recreate the figures in the article. Takes about 6 minutes. 

### sessionInfo.txt

The output of sessionInfo() when figurescript was run.

### overviewplot.r

Figure 1. in the manuscript. 

### reanalysis/  

The individual data sets that are used to illustrate the adverse effect of ComBat.
Consist of a figure-script that loads the data, process with or without ComBat, run a significance analysis and creates a plot of the p-values. Figures 2. and 3 in the manuscript.

### tex/

The tex files and related files needed to compile the manuscript as a pdf

### Additional resources

In our work, many more analyses than presented in the above figure-scripts were performed. These are more detailed and go beyond the scope of our article. But for the especially interested some of them are made available as rmarkdown reports here:
https://github.com/vegardny/batchnormwarning_reports.git
