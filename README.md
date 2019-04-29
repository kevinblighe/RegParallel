Standard regression functions in R enabled for parallel processing over large data-frames
================
Kevin Blighe
2019-04-29

-   [Introduction](#introduction)
-   [Installation](#installation)
    -   [1. Download the package from Bioconductor](#download-the-package-from-bioconductor)
    -   [2. Load the package into R session](#load-the-package-into-r-session)
-   [Quick start](#quick-start)
    -   [Perform the most basic logistic regression analysis](#perform-the-most-basic-logistic-regression-analysis)
    -   [Perform a basic linear regression](#perform-a-basic-linear-regression)
    -   [Perform the most basic negative binomial logistic regression analysis](#perform-the-most-basic-negative-binomial-logistic-regression-analysis)
    -   [Survival analysis via Cox Proportional Hazards regression](#survival-analysis-via-cox-proportional-hazards-regression)
    -   [Perform a conditional logistic regression](#perform-a-conditional-logistic-regression)
-   [Advanced features](#advanced-features)
    -   [Speed up processing](#speed-up-processing)
        -   [~2000 tests; blocksize, 500; cores, 2; nestedParallel, TRUE](#tests-blocksize-500-cores-2-nestedparallel-true)
        -   [~2000 tests; blocksize, 500; cores, 2; nestedParallel, FALSE](#tests-blocksize-500-cores-2-nestedparallel-false)
        -   [~40000 tests; blocksize, 2000; cores, 2; nestedParallel, TRUE](#tests-blocksize-2000-cores-2-nestedparallel-true)
        -   [~40000 tests; blocksize, 2000; cores, 2; nestedParallel, FALSE](#tests-blocksize-2000-cores-2-nestedparallel-false)
        -   [~40000 tests; blocksize, 5000; cores, 3; nestedParallel, TRUE](#tests-blocksize-5000-cores-3-nestedparallel-true)
    -   [Modify confidence intervals](#modify-confidence-intervals)
    -   [Remove some terms from output / include the intercept](#remove-some-terms-from-output-include-the-intercept)
-   [Acknowledgments](#acknowledgments)
-   [Session info](#session-info)
-   [References](#references)

Introduction
============

In many analyses, a large amount of variables have to be tested independently against the trait/endpoint of interest, and also adjusted for covariates and confounding factors at the same time. The major bottleneck in these is the amount of time that it takes to complete these analyses.

With <i>RegParallel</i>, a large number of tests can be performed simultaneously. On a 12-core system, 144 variables can be tested simultaneously, with 1000s of variables processed in a matter of seconds via 'nested' parallel processing.

Works for logistic regression, linear regression, conditional logistic regression, Cox proportional hazards and survival models, Bayesian logistic regression, and negative binomial regression.

Installation
============

1. Download the package from Bioconductor
-----------------------------------------

``` r
    if (!requireNamespace('BiocManager', quietly = TRUE))
        install.packages('BiocManager')
        BiocManager::install('RegParallel')
```

Note: to install development version:

``` r
    devtools::install_github('kevinblighe/RegParallel')
```

2. Load the package into R session
----------------------------------

``` r
    library(RegParallel)
```

Quick start
===========

For this quick start, we will follow the tutorial (from Section 3.1) of [RNA-seq workflow: gene-level exploratory analysis and differential expression](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html). Specifically, we will load the 'airway' data, where different airway smooth muscle cells were treated with dexamethasone.

``` r
    library(airway)
    library(magrittr)

    data('airway')
    airway$dex %<>% relevel('untrt')
```

Normalise the raw counts in <i>DESeq2</i> and produce regularised log counts:

``` r
    library(DESeq2)

    dds <- DESeqDataSet(airway, design = ~ dex + cell)
    dds <- DESeq(dds, betaPrior = FALSE)
    rlogcounts <- assay(rlog(dds, blind = FALSE))
    rlogdata <- data.frame(colData(airway), t(rlogcounts))
```

Perform the most basic logistic regression analysis
---------------------------------------------------

Here, we fit a binomial logistic regression model to the data via <i>glmParallel</i>, with dexamethasone as the dependent variable.

``` r
    res1 <- RegParallel(
      data = rlogdata[ ,1:3000],
      formula = 'dex ~ [*]',
      FUN = function(formula, data)
        glm(formula = formula,
          data = data,
          family = binomial(link = 'logit')),
      FUNtype = 'glm',
      variables = colnames(rlogdata)[10:3000])

    res1[order(res1$P, decreasing=FALSE),]
```

    ##              Variable            Term       Beta StandardError
    ##    1: ENSG00000095464 ENSG00000095464   43.27934  2.593463e+01
    ##    2: ENSG00000071859 ENSG00000071859   12.96251  7.890287e+00
    ##    3: ENSG00000069812 ENSG00000069812  -44.37139  2.704021e+01
    ##    4: ENSG00000072415 ENSG00000072415  -19.90841  1.227527e+01
    ##    5: ENSG00000073921 ENSG00000073921   14.59470  8.999831e+00
    ##   ---                                                         
    ## 2817: ENSG00000068831 ENSG00000068831  110.84893  2.729072e+05
    ## 2818: ENSG00000069020 ENSG00000069020 -186.45744  4.603615e+05
    ## 2819: ENSG00000083642 ENSG00000083642 -789.55666  1.951104e+06
    ## 2820: ENSG00000104331 ENSG00000104331  394.14700  9.749138e+05
    ## 2821: ENSG00000083097 ENSG00000083097 -217.48873  5.398191e+05
    ##                   Z          P            OR      ORlower      ORupper
    ##    1:  1.6687854476 0.09515991  6.251402e+18 5.252646e-04 7.440065e+40
    ##    2:  1.6428433092 0.10041536  4.261323e+05 8.190681e-02 2.217017e+12
    ##    3: -1.6409412536 0.10080961  5.367228e-20 5.165170e-43 5.577191e+03
    ##    4: -1.6218306224 0.10483962  2.258841e-09 8.038113e-20 6.347711e+01
    ##    5:  1.6216635641 0.10487540  2.179701e+06 4.761313e-02 9.978541e+13
    ##   ---                                                                 
    ## 2817:  0.0004061781 0.99967592  1.383811e+48 0.000000e+00           NA
    ## 2818: -0.0004050239 0.99967684  1.053326e-81 0.000000e+00           NA
    ## 2819: -0.0004046717 0.99967712  0.000000e+00 0.000000e+00           NA
    ## 2820:  0.0004042891 0.99967742 1.499223e+171 0.000000e+00           NA
    ## 2821: -0.0004028919 0.99967854  3.514358e-95 0.000000e+00           NA

Perform a basic linear regression
---------------------------------

Here, we will perform the linear regression using both <i>glmParallel</i> and <i>lmParallel</i>. We will appreciate that a linear regression is the same using either function with the default settings.

Regularised log counts from our <i>DESeq2</i> data will be used.

``` r
  rlogdata <- rlogdata[ ,1:2000]

  res2 <- RegParallel(
    data = rlogdata,
    formula = '[*] ~ cell',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = colnames(rlogdata)[10:ncol(rlogdata)],
    p.adjust = "none")

  res3 <- RegParallel(
    data = rlogdata,
    formula = '[*] ~ cell',
    FUN = function(formula, data)
      lm(formula = formula,
        data = data),
    FUNtype = 'lm',
    variables = colnames(rlogdata)[10:ncol(rlogdata)],
    p.adjust = "none")

  subset(res2, P<0.05)
```

    ##             Variable        Term        Beta StandardError          t
    ##   1: ENSG00000001461 cellN061011 -0.46859875    0.10526111  -4.451775
    ##   2: ENSG00000001461 cellN080611 -0.84020922    0.10526111  -7.982143
    ##   3: ENSG00000001461  cellN61311 -0.87778101    0.10526111  -8.339082
    ##   4: ENSG00000001561 cellN080611 -1.71802758    0.13649920 -12.586357
    ##   5: ENSG00000001561  cellN61311 -1.05328889    0.13649920  -7.716448
    ##  ---                                                                 
    ## 519: ENSG00000092108 cellN061011 -0.12721659    0.01564082  -8.133625
    ## 520: ENSG00000092108  cellN61311 -0.12451203    0.01564082  -7.960708
    ## 521: ENSG00000092148 cellN080611 -0.34988071    0.10313461  -3.392467
    ## 522: ENSG00000092200 cellN080611  0.05906656    0.01521063   3.883241
    ## 523: ENSG00000092208 cellN080611 -0.28587683    0.08506716  -3.360602
    ##                 P        OR   ORlower   ORupper
    ##   1: 0.0112313246 0.6258787 0.5092039 0.7692873
    ##   2: 0.0013351958 0.4316202 0.3511586 0.5305181
    ##   3: 0.0011301853 0.4157043 0.3382098 0.5109554
    ##   4: 0.0002293465 0.1794197 0.1373036 0.2344544
    ##   5: 0.0015182960 0.3487887 0.2669157 0.4557753
    ##  ---                                           
    ## 519: 0.0012429963 0.8805429 0.8539591 0.9079544
    ## 520: 0.0013489163 0.8829276 0.8562718 0.9104133
    ## 521: 0.0274674209 0.7047722 0.5757851 0.8626549
    ## 522: 0.0177922771 1.0608458 1.0296864 1.0929482
    ## 523: 0.0282890537 0.7513552 0.6359690 0.8876762

``` r
  subset(res3, P<0.05)
```

    ##             Variable        Term        Beta StandardError          t
    ##   1: ENSG00000001461 cellN061011 -0.46859875    0.10526111  -4.451775
    ##   2: ENSG00000001461 cellN080611 -0.84020922    0.10526111  -7.982143
    ##   3: ENSG00000001461  cellN61311 -0.87778101    0.10526111  -8.339082
    ##   4: ENSG00000001561 cellN080611 -1.71802758    0.13649920 -12.586357
    ##   5: ENSG00000001561  cellN61311 -1.05328889    0.13649920  -7.716448
    ##  ---                                                                 
    ## 519: ENSG00000092108 cellN061011 -0.12721659    0.01564082  -8.133625
    ## 520: ENSG00000092108  cellN61311 -0.12451203    0.01564082  -7.960708
    ## 521: ENSG00000092148 cellN080611 -0.34988071    0.10313461  -3.392467
    ## 522: ENSG00000092200 cellN080611  0.05906656    0.01521063   3.883241
    ## 523: ENSG00000092208 cellN080611 -0.28587683    0.08506716  -3.360602
    ##                 P        OR   ORlower   ORupper
    ##   1: 0.0112313246 0.6258787 0.5092039 0.7692873
    ##   2: 0.0013351958 0.4316202 0.3511586 0.5305181
    ##   3: 0.0011301853 0.4157043 0.3382098 0.5109554
    ##   4: 0.0002293465 0.1794197 0.1373036 0.2344544
    ##   5: 0.0015182960 0.3487887 0.2669157 0.4557753
    ##  ---                                           
    ## 519: 0.0012429963 0.8805429 0.8539591 0.9079544
    ## 520: 0.0013489163 0.8829276 0.8562718 0.9104133
    ## 521: 0.0274674209 0.7047722 0.5757851 0.8626549
    ## 522: 0.0177922771 1.0608458 1.0296864 1.0929482
    ## 523: 0.0282890537 0.7513552 0.6359690 0.8876762

Perform the most basic negative binomial logistic regression analysis
---------------------------------------------------------------------

Here, we will utilise normalised, unlogged counts from <i>DESeq2</i>. Unlogged counts in RNA-seq naturally follow a negative binomial / Poisson-like distribution. <i>glm.nbParallel</i> will be used.

``` r
    nbcounts <- round(counts(dds, normalized = TRUE), 0)
    nbdata <- data.frame(colData(airway), t(nbcounts))

    res4 <- RegParallel(
      data = nbdata[ ,1:3000],
      formula = '[*] ~ dex',
      FUN = function(formula, data)
        glm.nb(formula = formula,
          data = data),
      FUNtype = 'glm.nb',
      variables = colnames(nbdata)[10:3000],
    p.adjust = "fdr")

    res4[order(res4$Theta, decreasing = TRUE),]
```

    ##              Variable   Term         Beta StandardError          Z
    ##    1: ENSG00000102226 dextrt -0.139286172    0.01862107 -7.4800288
    ##    2: ENSG00000102910 dextrt -0.030278805    0.01582026 -1.9139254
    ##    3: ENSG00000063601 dextrt -0.002822867    0.02656550 -0.1062607
    ##    4: ENSG00000083642 dextrt -0.128217080    0.02412001 -5.3157976
    ##    5: ENSG00000023041 dextrt  0.134325359    0.02728332  4.9233512
    ##   ---                                                             
    ## 2817: ENSG00000029559 dextrt -1.386294361    1.73450917 -0.7992430
    ## 2818: ENSG00000006128 dextrt  1.386294361    1.73450917  0.7992430
    ## 2819: ENSG00000101197 dextrt -1.386294361    1.73450917 -0.7992430
    ## 2820: ENSG00000102109 dextrt -0.318453731    1.39587712 -0.2281388
    ## 2821: ENSG00000069122 dextrt  0.552068582    1.61318296  0.3422232
    ##                  P        Theta      SEtheta  2xLogLik Dispersion
    ##    1: 7.430632e-14 1.987085e+08 1.508873e+10 -73.91315          1
    ##    2: 5.562968e-02 1.124953e+08 6.192688e+09 -78.06984          1
    ##    3: 9.153755e-01 6.328822e+07 4.095363e+09 -68.80847          1
    ##    4: 1.061911e-07 3.111265e+07 1.578546e+09 -72.68055          1
    ##    5: 8.507456e-07 3.003015e+07 1.594172e+09 -70.09368          1
    ##   ---                                                            
    ## 2817: 4.241495e-01 2.843297e-01 3.514737e-01 -15.24150          1
    ## 2818: 4.241495e-01 2.843297e-01 3.514737e-01 -15.24150          1
    ## 2819: 4.241495e-01 2.843297e-01 3.514737e-01 -15.24150          1
    ## 2820: 8.195383e-01 2.664530e-01 1.500629e-01 -41.63690          1
    ## 2821: 7.321830e-01 1.984580e-01 1.233914e-01 -37.96063          1
    ##              OR    ORlower     ORupper     P.adjust
    ##    1: 0.8699790 0.83880014   0.9023169 2.624401e-12
    ##    2: 0.9701750 0.94055425   1.0007286 1.401173e-01
    ##    3: 0.9971811 0.94658900   1.0504772 9.941050e-01
    ##    4: 0.8796624 0.83904459   0.9222465 1.361660e-06
    ##    5: 1.1437649 1.08420938   1.2065918 8.333172e-06
    ##   ---                                              
    ## 2817: 0.2500000 0.00834686   7.4878457 6.161820e-01
    ## 2818: 4.0000000 0.13354976 119.8055313 6.161820e-01
    ## 2819: 0.2500000 0.00834686   7.4878457 6.161820e-01
    ## 2820: 0.7272727 0.04715465  11.2168280 9.303492e-01
    ## 2821: 1.7368421 0.07355573  41.0113595 8.703426e-01

Survival analysis via Cox Proportional Hazards regression
---------------------------------------------------------

For this example, we will load breast cancer gene expression data with recurrence free survival (RFS) from [Gene Expression Profiling in Breast Cancer: Understanding the Molecular Basis of Histologic Grade To Improve Prognosis](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2990). Specifically, we will encode each gene's expression into Low|Mid|High based on Z-scores and compare these against RFS while adjusting for tumour grade in a Cox Proportional Hazards model.

NB - the tutorial for this section can be found in my Biostars post: [Tutorial: Survival analysis with gene expression](https://www.biostars.org/p/344233/)

Perform a conditional logistic regression
-----------------------------------------

In this example, we will re-use the Cox data for the purpose of performing conditional logistic regression with tumour grade as our grouping / matching factor. For this example, we will use ER status as the dependent variable and also adjust for age.

``` r
  x <- exprs(gset[[1]])
  x <- x[-grep('^AFFX', rownames(x)),]
  x <- scale(x)
  x <- x[,which(colnames(x) %in% rownames(metadata))]

  coxdata <- data.frame(metadata, t(x))

  colnames(coxdata)[1:8] <- c('Age', 'Distant.RFS', 'ER',
    'GGI', 'Grade', 'Node',
    'Size', 'Time.RFS')

  coxdata$Age <- as.numeric(gsub('^KJ', '', coxdata$Age))
  coxdata$Grade <- factor(coxdata$Grade, levels = c(1, 2, 3))
  coxdata$ER <- as.numeric(coxdata$ER)
  coxdata <- coxdata[!is.na(coxdata$ER),]

  res6 <- RegParallel(
    data = coxdata,
    formula = 'ER ~ [*] + Age + strata(Grade)',
    FUN = function(formula, data)
      clogit(formula = formula,
        data = data,
        method = 'breslow'),
    FUNtype = 'clogit',
    variables = colnames(coxdata)[9:ncol(coxdata)],
    blocksize = 2000)

  subset(res6, P < 0.01)

  getBM(mart = mart,
    attributes = c('affy_hg_u133a',
      'ensembl_gene_id',
      'gene_biotype',
      'external_gene_name'),
    filter = 'affy_hg_u133a',
    values = c('204667_at',
      '205225_at',
      '207813_s_at',
      '212108_at',
      '219497_s_at'),
    uniqueRows=TRUE)
```

Oestrogen receptor (<i>ESR1</i>) comes out top - makes sense! Also, although 204667\_at is not listed in <i>biomaRt</i>, it overlaps an exon of <i>FOXA1</i>, which also makes sense in relation to oestrogen signalling.

Advanced features
=================

Advanced features include the ability to modify block size, choose different numbers of cores, enable 'nested' parallel processing, modify limits for confidence intervals, and exclude certain model terms from output.

Speed up processing
-------------------

First create some test data for the purpose of benchmarking:

``` r
  options(scipen=10)
  options(digits=6)

  # create a data-matrix of 20 x 60000 (rows x cols) random numbers
  col <- 60000
  row <- 20
  mat <- matrix(
    rexp(col*row, rate = .1),
    ncol = col)

  # add fake gene and sample names
  colnames(mat) <- paste0('gene', 1:ncol(mat))

  rownames(mat) <- paste0('sample', 1:nrow(mat))

  # add some fake metadata
  modelling <- data.frame(
    cell = rep(c('B', 'T'), nrow(mat) / 2),
    group = c(rep(c('treatment'), nrow(mat) / 2), rep(c('control'), nrow(mat) / 2)),
    dosage = t(data.frame(matrix(rexp(row, rate = 1), ncol = row))),
    mat,
    row.names = rownames(mat))
```

### ~2000 tests; blocksize, 500; cores, 2; nestedParallel, TRUE

With 2 cores instead of the default of 3, coupled with nestedParallel being enabled, a total of 2 x 2 = 4 threads will be used.

``` r
  df <- modelling[ ,1:2000]
  variables <- colnames(df)[4:ncol(df)]

  ptm <- proc.time()

  res <- RegParallel(
    data = df,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 500,
    cores = 2,
    nestedParallel = TRUE,
    p.adjust = "BY")

  proc.time() - ptm
```

    ##    user  system elapsed 
    ##  16.284   5.068  16.951

### ~2000 tests; blocksize, 500; cores, 2; nestedParallel, FALSE

``` r
  df <- modelling[ ,1:2000]
  variables <- colnames(df)[4:ncol(df)]

  ptm <- proc.time()

  res <- RegParallel(
    data = df,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 500,
    cores = 2,
    nestedParallel = FALSE,
    p.adjust = "BY")

  proc.time() - ptm
```

    ##    user  system elapsed 
    ##  28.672   5.900  14.585

Focusing on the elapsed time (as system time only reports time from the last core that finished), we can see that nested processing has negligible improvement or may actually be slower under certain conditions when tested over a small number of variables. This is likely due to the system being slowed by simply managing the larger number of threads. Nested processing's benefits can only be gained when processing a large number of variables:

### ~40000 tests; blocksize, 2000; cores, 2; nestedParallel, TRUE

``` r
  df <- modelling[ ,1:40000]
  variables <- colnames(df)[4:ncol(df)]

  system.time(RegParallel(
    data = df,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 2000,
    cores = 2,
    nestedParallel = TRUE))
```

    ##    user  system elapsed 
    ##  364.88   35.68  306.40

### ~40000 tests; blocksize, 2000; cores, 2; nestedParallel, FALSE

``` r
  df <- modelling[,1:40000]
  variables <- colnames(df)[4:ncol(df)]

  system.time(RegParallel(
    data = df,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 2000,
    cores = 2,
    nestedParallel = FALSE))
```

    ##    user  system elapsed 
    ##  320.68    1.54  331.22

Performance is system-dependent and even increasing cores may not result in huge gains in time. Performance is a trade-off between cores, forked threads, blocksize, and the number of terms in each model.

### ~40000 tests; blocksize, 5000; cores, 3; nestedParallel, TRUE

In this example, we choose a large blocksize and 3 cores. With nestedParallel enabled, this translates to 9 simultaneous threads.

``` r
  df <- modelling[,1:40000]
  variables <- colnames(df)[4:ncol(df)]

  system.time(RegParallel(
    data = df,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 5000,
    cores = 3,
    nestedParallel = TRUE))
```

    ##    user  system elapsed 
    ## 734.288  29.580 481.315

Modify confidence intervals
---------------------------

``` r
  df <- modelling[ ,1:500]
  variables <- colnames(df)[4:ncol(df)]

  # 99% confidfence intervals
  RegParallel(
    data = df,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 150,
    cores = 3,
    nestedParallel = TRUE,
    conflevel = 99)
```

    ##       Variable         Term       Beta StandardError         Z         P
    ##    1:    gene1        gene1  0.0369753     0.0838422  0.441010 0.6592057
    ##    2:    gene1 cellB:dosage -3.3954361     3.7191779 -0.912953 0.3612671
    ##    3:    gene1 cellT:dosage -0.7235396     0.9073731 -0.797400 0.4252186
    ##    4:    gene2        gene2 -0.0203082     0.0594907 -0.341367 0.7328273
    ##    5:    gene2 cellB:dosage -3.4132639     3.5460181 -0.962562 0.3357671
    ##   ---                                                                   
    ## 1487:  gene496 cellB:dosage -1.4091922     2.5665623 -0.549058 0.5829655
    ## 1488:  gene496 cellT:dosage -1.6507933     1.1916116 -1.385345 0.1659470
    ## 1489:  gene497      gene497  0.1779883     0.1029215  1.729359 0.0837448
    ## 1490:  gene497 cellB:dosage -6.9433365     4.8094909 -1.443674 0.1488307
    ## 1491:  gene497 cellT:dosage -3.0020996     1.8662221 -1.608651 0.1076927
    ##                OR          ORlower   ORupper
    ##    1: 1.037667345 0.83611597240022   1.28780
    ##    2: 0.033525930 0.00000231661476 485.18554
    ##    3: 0.485032411 0.04685124609439   5.02135
    ##    4: 0.979896646 0.84067834779231   1.14217
    ##    5: 0.032933532 0.00000355483746 305.11030
    ##   ---                                       
    ## 1487: 0.244340580 0.00032874813411 181.60504
    ## 1488: 0.191897615 0.00891356072518   4.13131
    ## 1489: 1.194811336 0.91656711941876   1.55752
    ## 1490: 0.000965044 0.00000000402088 231.61838
    ## 1491: 0.049682646 0.00040599526609   6.07979

``` r
  # 95% confidfence intervals (default)
  RegParallel(
    data = df,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 150,
    cores = 3,
    nestedParallel = TRUE,
    conflevel = 95)
```

    ##       Variable         Term       Beta StandardError         Z         P
    ##    1:    gene1        gene1  0.0369753     0.0838422  0.441010 0.6592057
    ##    2:    gene1 cellB:dosage -3.3954361     3.7191779 -0.912953 0.3612671
    ##    3:    gene1 cellT:dosage -0.7235396     0.9073731 -0.797400 0.4252186
    ##    4:    gene2        gene2 -0.0203082     0.0594907 -0.341367 0.7328273
    ##    5:    gene2 cellB:dosage -3.4132639     3.5460181 -0.962562 0.3357671
    ##   ---                                                                   
    ## 1487:  gene496 cellB:dosage -1.4091922     2.5665623 -0.549058 0.5829655
    ## 1488:  gene496 cellT:dosage -1.6507933     1.1916116 -1.385345 0.1659470
    ## 1489:  gene497      gene497  0.1779883     0.1029215  1.729359 0.0837448
    ## 1490:  gene497 cellB:dosage -6.9433365     4.8094909 -1.443674 0.1488307
    ## 1491:  gene497 cellT:dosage -3.0020996     1.8662221 -1.608651 0.1076927
    ##                OR         ORlower  ORupper
    ##    1: 1.037667345 0.8804233162184  1.22300
    ##    2: 0.033525930 0.0000228881581 49.10784
    ##    3: 0.485032411 0.0819244272252  2.87163
    ##    4: 0.979896646 0.8720505599278  1.10108
    ##    5: 0.032933532 0.0000315691033 34.35693
    ##   ---                                     
    ## 1487: 0.244340580 0.0015971061783 37.38156
    ## 1488: 0.191897615 0.0185681189796  1.98322
    ## 1489: 1.194811336 0.9765452572486  1.46186
    ## 1490: 0.000965044 0.0000000777501 11.97825
    ## 1491: 0.049682646 0.0012813672456  1.92635

Remove some terms from output / include the intercept
-----------------------------------------------------

``` r
  # remove terms but keep Intercept
  RegParallel(
    data = df,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 150,
    cores = 3,
    nestedParallel = TRUE,
    conflevel = 95,
    excludeTerms = c('cell', 'dosage'),
    excludeIntercept = FALSE)
```

    ##      Variable        Term       Beta StandardError          Z         P
    ##   1:    gene1 (Intercept)  0.4668797     1.0671820  0.4374884 0.6617572
    ##   2:    gene1       gene1  0.0369753     0.0838422  0.4410101 0.6592057
    ##   3:    gene2 (Intercept)  0.9113724     0.8713520  1.0459291 0.2955938
    ##   4:    gene2       gene2 -0.0203082     0.0594907 -0.3413670 0.7328273
    ##   5:    gene3 (Intercept)  0.6570693     0.8197455  0.8015527 0.4228117
    ##  ---                                                                   
    ## 990:  gene495     gene495  0.0466219     0.0420685  1.1082390 0.2677586
    ## 991:  gene496 (Intercept) -0.8248775     1.0523229 -0.7838635 0.4331202
    ## 992:  gene496     gene496  0.1774098     0.1106534  1.6032926 0.1088701
    ## 993:  gene497 (Intercept) -0.0647892     0.9899981 -0.0654437 0.9478207
    ## 994:  gene497     gene497  0.1779883     0.1029215  1.7293593 0.0837448
    ##            OR   ORlower  ORupper
    ##   1: 1.595010 0.1969592 12.91666
    ##   2: 1.037667 0.8804233  1.22300
    ##   3: 2.487734 0.4509287 13.72462
    ##   4: 0.979897 0.8720506  1.10108
    ##   5: 1.929130 0.3868948  9.61901
    ##  ---                            
    ## 990: 1.047726 0.9648036  1.13778
    ## 991: 0.438289 0.0557213  3.44746
    ## 992: 1.194120 0.9613018  1.48333
    ## 993: 0.937265 0.1346401  6.52455
    ## 994: 1.194811 0.9765453  1.46186

``` r
  # remove everything but the variable being tested
  RegParallel(
    data = df,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 150,
    cores = 3,
    nestedParallel = TRUE,
    conflevel = 95,
    excludeTerms = c('cell', 'dosage'),
    excludeIntercept = TRUE)
```

    ##      Variable    Term       Beta StandardError         Z         P
    ##   1:    gene1   gene1  0.0369753     0.0838422  0.441010 0.6592057
    ##   2:    gene2   gene2 -0.0203082     0.0594907 -0.341367 0.7328273
    ##   3:    gene3   gene3  0.0412457     0.0780303  0.528586 0.5970925
    ##   4:    gene4   gene4  0.0257413     0.0507055  0.507664 0.6116890
    ##   5:    gene5   gene5  0.1604644     0.0787428  2.037828 0.0415671
    ##  ---                                                              
    ## 493:  gene493 gene493 -0.3535213     0.2660887 -1.328584 0.1839851
    ## 494:  gene494 gene494  0.0245676     0.0931056  0.263869 0.7918811
    ## 495:  gene495 gene495  0.0466219     0.0420685  1.108239 0.2677586
    ## 496:  gene496 gene496  0.1774098     0.1106534  1.603293 0.1088701
    ## 497:  gene497 gene497  0.1779883     0.1029215  1.729359 0.0837448
    ##            OR  ORlower ORupper
    ##   1: 1.037667 0.880423 1.22300
    ##   2: 0.979897 0.872051 1.10108
    ##   3: 1.042108 0.894321 1.21432
    ##   4: 1.026076 0.929006 1.13329
    ##   5: 1.174056 1.006150 1.36998
    ##  ---                          
    ## 493: 0.702211 0.416843 1.18294
    ## 494: 1.024872 0.853922 1.23005
    ## 495: 1.047726 0.964804 1.13778
    ## 496: 1.194120 0.961302 1.48333
    ## 497: 1.194811 0.976545 1.46186

Acknowledgments
===============

*RegParallel* would not exist were it not for initial contributions from:

[Jessica Lasky-Su](https://connects.catalyst.harvard.edu/Profiles/display/Person/79780), [Myles Lewis](https://www.qmul.ac.uk/whri/people/academic-staff/items/lewismyles.html), [Michael Barnes](https://www.qmul.ac.uk/whri/people/academic-staff/items/barnesmichael.html)

Thanks also to Horacio Montenegro and GenoMax for testing.

Session info
============

``` r
sessionInfo()
```

    ## R version 3.6.0 (2019-04-26)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/atlas-base/atlas/libblas.so.3.0
    ## LAPACK: /usr/lib/atlas-base/atlas/liblapack.so.3.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=pt_BR.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=pt_BR.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=pt_BR.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] DESeq2_1.23.10              magrittr_1.5               
    ##  [3] airway_1.3.0                SummarizedExperiment_1.13.0
    ##  [5] DelayedArray_0.9.9          BiocParallel_1.17.19       
    ##  [7] matrixStats_0.54.0          Biobase_2.43.1             
    ##  [9] GenomicRanges_1.35.1        GenomeInfoDb_1.19.3        
    ## [11] IRanges_2.17.5              S4Vectors_0.21.24          
    ## [13] BiocGenerics_0.29.2         RegParallel_1.1.3          
    ## [15] arm_1.10-1                  lme4_1.1-21                
    ## [17] Matrix_1.2-17               MASS_7.3-51.4              
    ## [19] survival_2.44-1.1           stringr_1.4.0              
    ## [21] data.table_1.12.2           doParallel_1.0.14          
    ## [23] iterators_1.0.10            foreach_1.4.4              
    ## [25] knitr_1.22                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit64_0.9-7            splines_3.6.0          Formula_1.2-3         
    ##  [4] assertthat_0.2.1       latticeExtra_0.6-28    blob_1.1.1            
    ##  [7] GenomeInfoDbData_1.2.1 yaml_2.2.0             RSQLite_2.1.1         
    ## [10] backports_1.1.4        pillar_1.3.1           lattice_0.20-38       
    ## [13] glue_1.3.1             digest_0.6.18          RColorBrewer_1.1-2    
    ## [16] XVector_0.23.2         checkmate_1.9.1        minqa_1.2.4           
    ## [19] colorspace_1.4-1       htmltools_0.3.6        plyr_1.8.4            
    ## [22] XML_3.98-1.19          pkgconfig_2.0.2        genefilter_1.65.0     
    ## [25] zlibbioc_1.29.0        xtable_1.8-4           purrr_0.3.2           
    ## [28] scales_1.0.0           annotate_1.61.1        tibble_2.1.1          
    ## [31] htmlTable_1.13.1       ggplot2_3.1.1          nnet_7.3-12           
    ## [34] lazyeval_0.2.2         crayon_1.3.4           memoise_1.1.0         
    ## [37] evaluate_0.13          nlme_3.1-139           foreign_0.8-71        
    ## [40] tools_3.6.0            locfit_1.5-9.1         munsell_0.5.0         
    ## [43] cluster_2.0.8          AnnotationDbi_1.45.1   compiler_3.6.0        
    ## [46] rlang_0.3.4            grid_3.6.0             RCurl_1.95-4.12       
    ## [49] nloptr_1.2.1           rstudioapi_0.10        htmlwidgets_1.3       
    ## [52] bitops_1.0-6           base64enc_0.1-3        rmarkdown_1.12        
    ## [55] boot_1.3-22            gtable_0.3.0           codetools_0.2-16      
    ## [58] DBI_1.0.0              abind_1.4-5            R6_2.4.0              
    ## [61] gridExtra_2.3          dplyr_0.8.0.1          bit_1.1-14            
    ## [64] Hmisc_4.2-0            stringi_1.4.3          Rcpp_1.0.1            
    ## [67] geneplotter_1.61.0     rpart_4.1-15           acepack_1.4.1         
    ## [70] tidyselect_0.2.5       xfun_0.6               coda_0.19-2

References
==========

Blighe (2018)

Blighe, Kevin. 2018. “RegParallel: Standard regression functions in R enabled for parallel processing over large data-frames.” <https://github.com/kevinblighe>.
