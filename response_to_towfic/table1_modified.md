# Results

Results are produced by running table1_modified_code.R on the data set GEO:GSE61901.

## Modified Table 1

                                                                    |limma |anova |mixed |aov+err
------------------------------------------------------------------- | ---- | ---- | ---- | ------
1. Towfic et al.: technical replicates, CHIP.blocking               | 1474 |   NA | 1968 |   1749
2. Equivalent but simpler model                                     |   NA |   NA |   NA |   1749
3. Corrected models: technical replicates, CHIP.blocking            |   NA |   NA |    1 |      1
4. Same model, different implementation                             |   NA |   NA |    0 |     NA
5. One-way: replicates averaged, ComBat with treatment as covariate |  974 |  742 |   NA |     NA
6. Two-way: replicates averaged, CHIP blocking                      |   87 |    1 |   NA |     NA
7. Two-way: replicates averaged, ComBat normalized, CHIP blocking   |   86 |    1 |   NA |     NA
8. Result from Nygaard et al. article (GEO:GSEGSE40566)             |   11 |   NA |   NA |     NA

Consensus correlation from duplicateCorrelation: 0.09611996

## Comments on results and explanations:

**1. These were the original "alternative methods" presented in the letter by Towfic et al.**

limma: Linear model, lmFit, using consensus correlation from duplicateCorrelation.

mixed: Mixed random effects model using lmer from the lme4 package. Model comparison using likelihood ratio test, which is not appropriate as it does not correct for the effect of the random effects on the estimates. Uses ML as REML does not make sense in such a comparison, but which tend to give biased variance estimates.

aov+err: Uses aov with random sample effect misspecified as Error(1/sampid). See 2.

**2. Demonstration that misspecified model is equivalent to model without random effects term**

aov+err: Same model as in 1, but without random effects included. Results exactly the same.

**3. Corrected mixed random effects models**

mixed: Fits model with lmer using REML to get unbiased variance estimates. Uses lmerTest to estimate degrees of freedom: done using Satterthwaite's method.

aov+err: Random effects term specified as Error(sampid).

**4. Alternative mixed random effects model**

mixed: Uses lme from the nlme package. Failed on some probes returning NA instead.

**5. One-way ANOVA/LIMMA after ComBat batch adjustment as done by Towfic et al.**

This approach was demonstrated in our article to be unreliable as batch adjustments are performed without taking into regard the estimation error of the batch effects. 

**6. Two-way ANOVA/LIMMA blocking by CHIP (batch)**

This is the same analysis as in the letter by Towfic et al., but we have included the actual LIMMA analysis on GEO:GSE61901 rather than citing the result from our article which was on GEO:GSE40566.

**7. Same as 6 but applied to ComBat normalised data**

ComBat may also rescale data per batch, and results may thus differ.

**8. Result from our article where we analysed GEO:GSE40566 which did not have replicates included**

The was based on the same method as in 7.


## P values for first ten probes

**Original analyses (and equivalent model) which we consider incorrect or unreliable**

             |limma.dupcor |ranova.orig |ranova.equiv | lme4.orig
-------------|-------------|------------|-------------|----------
ILMN_1238828 |   0.7439403 |  0.7094641 |   0.7094641 | 0.6759463
ILMN_1214761 |   0.4716057 |  0.4645878 |   0.4645878 | 0.4119452
ILMN_2767246 |   0.2466421 |  0.2146228 |   0.2146228 | 0.1638481
ILMN_2844032 |   0.5653758 |  0.5483816 |   0.5483816 | 0.5007167
ILMN_2741117 |   0.5572201 |  0.5512430 |   0.5512430 | 0.5288311
ILMN_2670150 |   0.0000007 |  0.0000051 |   0.0000051 | 0.0002116
ILMN_1221475 |   0.3128733 |  0.2529046 |   0.2529046 | 0.1995927
ILMN_1234913 |   0.3761661 |  0.3017442 |   0.3017442 | 0.2466173
ILMN_2610645 |   0.7176044 |  0.6944818 |   0.6944818 | 0.6594528
ILMN_2883907 |   0.0235338 |  0.0146523 |   0.0146523 | 0.0152355
_*) Note that ranova.orig = ranova.equiv_

**One-way models after ComBat batch adjustments, which we have demonstrated to be an unreliable approach**

             | limma.oneway| anova.oneway
-------------|-------------|-------------
ILMN_1238828 |   0.9207860 |    0.9170674 
ILMN_1214761 |   0.2961804 |    0.3096932 
ILMN_2767246 |   0.3145038 |    0.3055694 
ILMN_2844032 |   0.6094450 |    0.5982844 
ILMN_2741117 |   0.6790801 |    0.6935695 
ILMN_2670150 |   0.0000540 |    0.0002452 
ILMN_1221475 |   0.5071345 |    0.4869344 
ILMN_1234913 |   0.3804040 |    0.3519286 
ILMN_2610645 |   0.6564963 |    0.6448149 
ILMN_2883907 |   0.2867952 |    0.3031055 

**Random effects models and two-way ANOVA**

             |    ranova |     lme4  |      nlme |anova.twoway
-------------|-----------|-----------|-----------|------------
ILMN_1238828 | 0.7056701 | 0.7094641 | 0.7112863 |   0.7056701
ILMN_1214761 | 0.5066785 | 0.5066785 | 0.5066785 |   0.5066785
ILMN_2767246 | 0.1102511 | 0.2146228 | 0.2213085 |   0.1102511
ILMN_2844032 | 0.5096139 | 0.5483816 | 0.5515190 |   0.5096139
ILMN_2741117 | 0.6288526 | 0.6288526 | 0.6288526 |   0.6288526
ILMN_2670150 | 0.0044562 | 0.0044562 | 0.0044562 |   0.0044562
ILMN_1221475 | 0.2669480 | 0.2669480 | 0.2669480 |   0.2669480
ILMN_1234913 | 0.3350278 | 0.3350278 | 0.3350278 |   0.3350278
ILMN_2610645 | 0.6833171 | 0.6944818 | 0.6964139 |   0.6833171
ILMN_2883907 | 0.0625173 | 0.0625173 | 0.0625173 |   0.0625173
_*) Note that ranova = anova.twoway = lme4 = nlme in most rows, but with some rows being exceptions (still with ranova = anova.twoway) which we suspect is due to implementation differences (eg negative variance estimates)_

**Remaining two-way models**

             | limma.twoway | limma.combat |anova.combat
-------------|--------------|--------------|-----------
ILMN_1238828 |   0.7343685  |   0.6168479  |  0.5836370
ILMN_1214761 |   0.4811775  |   0.3154480  |  0.3428109
ILMN_2767246 |   0.1481677  |   0.3035237  |  0.2874723
ILMN_2844032 |   0.5137711  |   0.7093633  |  0.7010447
ILMN_2741117 |   0.6035043  |   0.4754986  |  0.5095846
ILMN_2670150 |   0.0009892  |   0.0006757  |  0.0039771
ILMN_1221475 |   0.3071734  |   0.3102254  |  0.2698768
ILMN_1234913 |   0.3850121  |   0.5598159  |  0.5176542
ILMN_2610645 |   0.6983242  |   0.7515402  |  0.7362996
ILMN_2883907 |   0.0555117  |   0.0744223  |  0.0832484
