Supplementary information in response to the Letter to the Editor by Towfic et al. in Biostatistics
----------

In our article, we reanalyzed parts of a study by [Towfic et al.](https://dx.doi.org/10.1371/journal.pone.0083757) in order to illustrate our main point, that batch effect adjustment may lead to unreliable results. Towfic et al. maintain that their results are reliable, and wrote a [Letter to the Editor](https://dx.doi.org/10.1093/biostatistics/kxw031) where they present three alternative ways to analyze the data, claiming to obtain results consistent with their original findings. The three new alternative ways to analyze the data all handle batch as a fixed effect in the models,<i>i.e.</i> in a way we consider to be reliable. However, for other reasons, the implementation of the analysis is problematic for the first one, using duplicateCorrelation, and just wrong for the two other. In order to realize this, we had to dive deep into their code. And for others to also realize the errors, there is no way around getting into the details.


### Problem with alternative [analysis 1](https://github.com/immuneering/biostat_reply/blob/master/table1_code.R#L96): duplicateCorrelation+LIMMA

The first alternative method uses the function duplicateCorrelation to estimate the correlation
between the two strips, and then LIMMA's lmFit with correlation specified to adjust for the
correlation between them. This method applies a single correlation across all probes, not a
correlation fitted per probe. When the correlation between the two strips is fairly consistent across
probes, or even when it displays some degree of variation, this approach may be acceptable.

In this case, the consensus correlation is 0.096, which in effect means the two strips are considered to
be almost independent. The reason is that the majority of transcripts are unexpressed, and the probes for the unexpressed genes primarily
capture technical variation. However, for well-expressed transcripts, the consensus correlation
rises above 0.5, and yet higher for transcripts that display substantial variation between samples.

By applying a consensus correlation close to zero, LIMMA's lmFIT will treat the two strips as near independent. This makes
the analysis in effect similar to a two-way ANOVA, with the two strips included as independent observations. For unexpressed
transcripts, which are uncorrelated, this is ok. However, for expressed transcripts where the two strips are strongly correlated,
this is highly inappropriate.

In a more recent, albeit brief [forum post](https://support.bioconductor.org/p/87680) regarding an experiment similar in design to the one discussed here, Gordon Smyth (author of duplicateCorrelation) replies that the method will lead to too small p-values.


### Problem with alternative [analysis 2](https://github.com/immuneering/biostat_reply/blob/master/table1_code.R#L47): aov with random effects

The second alternative method, presented as a repeated measures model, is a mixed random
effects model using R's aov function. Unfortunately, the random sample effect is misspecified
as <i>Error(1/sampid)</i> ([line 48 of their script](https://github.com/immuneering/biostat_reply/blob/master/table1_code.R#L48), which has no effect on the model, leaving the two strips to be treated as
independent observations. The correct syntax for an error term per sample in aov would be
<i>Error(sampid)</i>, which makes the model equivalent to the two-way ANOVA on data where the
two strips have been averaged. Thus, instead of the 1749 differentially expressed probesets, a
correct analysis only finds 1.

R´s formulae syntax is admittedly intricate. One quick way to see that <i>Error(1/sampid)</i> is nonsensical is to replace it with <i>Error(1)</i> and rerun. Leaving <i>sampid</i> out will produce identical results, demonstrating that the technical replicate information, identified by <i>sampid</i>, is not taken into account as part of the differential expression analysis. Leaving out the whole Error term will yield the same results, but more code needs to be modified, since the results from aov are differently formatted.


### Problem with alternative [analysis 3](https://github.com/immuneering/biostat_reply/blob/master/table1_code.R#L65): lmer mixed random effects model

The third alternative method is a mixed random effects model using R's lmer function from the lme4 package,
fitting maximum likelihood models with and without the treatment effect, and comparing their fit using
a likelihood ratio test.

One factor critical for this particular case is that the maximum likelihood estimates will tend to
be biased: particularly for small samples or models with many degrees of freedom.
Another is that the method ends up using inappropriate degrees of freedom for the chi-square statistics:
<i>e.g.</i> see [Halekoh and Højsgaard (2014)](https://www.jstatsoft.org/htaccess.php?volume=59&type=i&issue=09&paper=true)
for a discussion. This can easily be verified by running the analyses using random data with no group effect: see
[random.R](random.R).

Instead, the lmerTest library may be used to estimate appropriate degrees of freedom
required for significance testing. Alternatively,
the function lme from the nlme may be used. Again, since there are exactly two replicates per
samples, this becomes essentially identical to the two-way ANOVA after averaging the strips,
providing 1 differentially expressed probe, instead of the 1968 they report: the technical difference is due
to lmer enforcing non-negative variance estimates.


### Common problem

We also note that the reanalyses by Towfic et al. in the letter, unlike their and our original analyses, filter out samples from the original study that are not part of the comparison of Copaxone versus the generic drug. This renders many batches with samples from only one treatment group, and more than half the samples not contributing to the estimation of the difference. Thus, instead of comparing 34 against 11 samples from 17 batches, the comparison is effectively of 12 against 9 samples from 7 batches.


### Corrected scripts and results

**[table1_modified_code.R](table1_modified_code.R):** R script, modified from the one used by Towfic et al., with corrected models and alternative models

**[table1_modified.md](table1_modified.md):** Results from R script with an extended version of Table 1 from Towfic et al., extended comments, and comparison of results from different models.


### Other

**[Illumina_correspondence.txt](illumina_correspondence.txt):** Copy of email correspondence with Illumina tech support regarding processing of strips from the MouseWG-6 v2.0 bead array chip.


### References

Towfics et al.´s [**original research article**](https://dx.doi.org/10.1371/journal.pone.0083757)  
Towfic, F., Funt, J. M., Fowler, K. D., Bakshi, S., Blaugrund, E., Artyomov, M. N.,
Hayden, M. R., Ladkani, D., Schwartz, R. and Zeskind, B. (2014, jan). Comparing the
biological impact of glatiramer acetate with the biological impact of a generic. PloS one 9(1),
e83757.
<br/>

[**Our article, with the reanlysis**](https://dx.doi.org/10.1093/biostatistics/kxv027)  
Nygaard V, Rødland EA, Hovig E. Methods that remove batch effects while
retaining group differences may lead to exaggerated confidence in downstream
analyses. Biostatistics. 2016 Jan;17(1):29-39. doi: 10.1093/biostatistics/kxv027.
PubMed PMID: 26272994; PubMed Central PMCID: PMC4679072.

<br/>
Towfics et al.´s [**Letter to the Editor**](https://dx.doi.org/10.1093/biostatistics/kxw031) and the [r-script](https://github.com/immuneering/biostat_reply) to reproduce the results presented.

<br/>
[**GEO:GSE40566**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40566)  
Data set cited (but not used) in the original article by Towfic et al. Also used in our reanalysis. Data from the two strips are combined, thus each sample has one measurement per probe.
<br/>

[**GEO:GSE61901**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61901)  
Data set used in the original article by Towfic et al. and in their Letter to the Editor. Data from the two strips are kept as separate observations to be used as technical replicates, or erroneously as biological replicates.

<br/>
<br/>
