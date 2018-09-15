# StatParallel
Standard regression, correlation, and other functions in R enabled for parallel processing over large data-frames.

In many analyses, a large amount of variables have to be tested independently against the trait/endpoint of interest, and also adjusted for covariates and confounding factors at the same time. The major botteleneck in these is the amount of time that it takes to complete these analyses.

With <i>StatParallel</i>, any number of tests can be performed simultaneously.  On a 12-core system, 144 variables can be tested simultaneously, with 1000s of variables processed in a matter of seconds.

Works for logistic regression, linear regression, conditional logistic regression, Cox proportional hazards models, ANOVA, and correlations. Also works for GWAS studies loaded into R as snpMatrix objects.


<h1>Credits</h1>
<ul>
  <li>Kevin Blighe (Brigham & Women's Hospital / Harvard Medical School)</li>
  <li>Jessica Lasky-Su (Brigham & Women's Hospital / Harvard Medical School)</li>
  <li>Juan Celedón (Children's Hospital of Pittsburgh)</li>
  <li>Qi Yan (Children's Hospital of Pittsburgh)</li>
