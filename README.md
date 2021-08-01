Standard regression functions in R enabled for parallel processing over
large data-frames
================
Kevin Blighe, Jessica Lasky-Su
2021-08-01

# Introduction

In many analyses, a large amount of variables have to be tested
independently against the trait/endpoint of interest, and also adjusted
for covariates and confounding factors at the same time. The major
bottleneck in these is the amount of time that it takes to complete
these analyses.

With <i>RegParallel</i>, a large number of tests can be performed
simultaneously. On a 12-core system, 144 variables can be tested
simultaneously, with 1000s of variables processed in a matter of seconds
via ‘nested’ parallel processing.

Works for logistic regression, linear regression, conditional logistic
regression, Cox proportional hazards and survival models, and Bayesian
logistic regression. Also caters for generalised linear models that
utilise survey weights created by the ‘survey’ CRAN package and that
utilise ‘survey::svyglm’.

# Installation

## 1\. Download the package from Bioconductor

``` r
    if (!requireNamespace('BiocManager', quietly = TRUE))
        install.packages('BiocManager')

    BiocManager::install('RegParallel')
```

Note: to install development version:

``` r
    devtools::install_github('kevinblighe/RegParallel')
```

## 2\. Load the package into R session

``` r
    library(RegParallel)
```

# Quick start

For this quick start, we will follow the tutorial (from Section 3.1) of
[RNA-seq workflow: gene-level exploratory analysis and differential
expression](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html).
Specifically, we will load the ‘airway’ data, where different airway
smooth muscle cells were treated with dexamethasone.

``` r
    library(airway)
    library(magrittr)

    data('airway')
    airway$dex %<>% relevel('untrt')
```

Normalise the raw counts in <i>DESeq2</i> and produce regularised log
expression levels:

``` r
    library(DESeq2)

    dds <- DESeqDataSet(airway, design = ~ dex + cell)
    dds <- DESeq(dds, betaPrior = FALSE)
    rldexpr <- assay(rlog(dds, blind = FALSE))
    rlddata <- data.frame(colData(airway), t(rldexpr))
```

## Perform the most basic logistic regression analysis

Here, we fit a binomial logistic regression model to the data via
<i>glmParallel</i>, with dexamethasone as the dependent variable.

``` r
    ## NOT RUN

    res1 <- RegParallel(
      data = rlddata[ ,1:3000],
      formula = 'dex ~ [*]',
      FUN = function(formula, data)
        glm(formula = formula,
          data = data,
          family = binomial(link = 'logit')),
      FUNtype = 'glm',
      variables = colnames(rlddata)[10:3000])

    res1[order(res1$P, decreasing=FALSE),]
```

## Perform a basic linear regression

Here, we will perform the linear regression using both
<i>glmParallel</i> and <i>lmParallel</i>. We will appreciate that a
linear regression is the same using either function with the default
settings.

Regularised log expression levels from our <i>DESeq2</i> data will be
used.

``` r
  rlddata <- rlddata[ ,1:2000]

  res2 <- RegParallel(
    data = rlddata,
    formula = '[*] ~ cell',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = colnames(rlddata)[10:ncol(rlddata)],
    p.adjust = "none")

  res3 <- RegParallel(
    data = rlddata,
    formula = '[*] ~ cell',
    FUN = function(formula, data)
      lm(formula = formula,
        data = data),
    FUNtype = 'lm',
    variables = colnames(rlddata)[10:ncol(rlddata)],
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

## Survival analysis via Cox Proportional Hazards regression

For this example, we will load breast cancer gene expression data with
recurrence free survival (RFS) from [Gene Expression Profiling in Breast
Cancer: Understanding the Molecular Basis of Histologic Grade To Improve
Prognosis](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2990).
Specifically, we will encode each gene’s expression into Low|Mid|High
based on Z-scores and compare these against RFS while adjusting for
tumour grade in a Cox Proportional Hazards model.

First, let’s read in and prepare the data:

``` r
  library(Biobase)
  library(GEOquery)
  # load series and platform data from GEO
  gset <- getGEO('GSE2990', GSEMatrix =TRUE, getGPL=FALSE)
  x <- exprs(gset[[1]])
  # remove Affymetrix control probes
  x <- x[-grep('^AFFX', rownames(x)),]
  # transform the expression data to Z scores
  x <- t(scale(t(x)))
  # extract information of interest from the phenotype data (pdata)
  idx <- which(colnames(pData(gset[[1]])) %in%
    c('age:ch1', 'distant rfs:ch1', 'er:ch1',
      'ggi:ch1', 'grade:ch1', 'node:ch1',
      'size:ch1', 'time rfs:ch1'))
  metadata <- data.frame(pData(gset[[1]])[,idx],
    row.names = rownames(pData(gset[[1]])))
  # remove samples from the pdata that have any NA value
  discard <- apply(metadata, 1, function(x) any(is.na(x)))
  metadata <- metadata[!discard,]
  # filter the Z-scores expression data to match the samples in our pdata
  x <- x[,which(colnames(x) %in% rownames(metadata))]
  # check that sample names match exactly between pdata and Z-scores 
  all((colnames(x) == rownames(metadata)) == TRUE)
```

    ## [1] TRUE

``` r
  # create a merged pdata and Z-scores object
  coxdata <- data.frame(metadata, t(x))
  # tidy column names
  colnames(coxdata)[1:8] <- c('Age', 'Distant.RFS', 'ER',
    'GGI', 'Grade', 'Node',
    'Size', 'Time.RFS')
  # prepare certain phenotypes
  coxdata$Age <- as.numeric(gsub('^KJ', '', coxdata$Age))
  coxdata$Distant.RFS <- as.numeric(coxdata$Distant.RFS)
  coxdata$ER <- factor(coxdata$ER, levels = c(0, 1))
  coxdata$Grade <- factor(coxdata$Grade, levels = c(1, 2, 3))
  coxdata$Time.RFS <- as.numeric(gsub('^KJX|^KJ', '', coxdata$Time.RFS))
```

With the data prepared, we can now apply a Cox Proportional Hazards
model independently for each probe in the dataset against RFS.

In this we also increase the default blocksize to 2000 in order to speed
up the analysis.

``` r
  library(survival)
  res5 <- RegParallel(
    data = coxdata,
    formula = 'Surv(Time.RFS, Distant.RFS) ~ [*]',
    FUN = function(formula, data)
      coxph(formula = formula,
        data = data,
        ties = 'breslow',
        singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = colnames(coxdata)[9:ncol(coxdata)],
    blocksize = 2000,
    p.adjust = "BH")
  res5 <- res5[!is.na(res5$P),]
  res5
```

    ##           Variable        Term          Beta StandardError             Z
    ##     1:  X1007_s_at  X1007_s_at  0.3780639987     0.3535022  1.0694811914
    ##     2:    X1053_at    X1053_at  0.1177398813     0.2275041  0.5175285346
    ##     3:     X117_at     X117_at  0.6265036787     0.6763106  0.9263549892
    ##     4:     X121_at     X121_at -0.6138126274     0.6166626 -0.9953783151
    ##     5:  X1255_g_at  X1255_g_at -0.2043297829     0.3983930 -0.5128849375
    ##    ---                                                                  
    ## 22211:   X91703_at   X91703_at -0.4124539527     0.4883759 -0.8445419981
    ## 22212: X91816_f_at X91816_f_at  0.0482030943     0.3899180  0.1236236554
    ## 22213:   X91826_at   X91826_at  0.0546751431     0.3319572  0.1647053850
    ## 22214:   X91920_at   X91920_at -0.6452125945     0.8534623 -0.7559942684
    ## 22215:   X91952_at   X91952_at -0.0001396044     0.7377681 -0.0001892254
    ##                P       LRT      Wald   LogRank        HR    HRlower  HRupper
    ##     1: 0.2848529 0.2826716 0.2848529 0.2848400 1.4594563 0.72994385 2.918050
    ##     2: 0.6047873 0.6085603 0.6047873 0.6046839 1.1249515 0.72024775 1.757056
    ##     3: 0.3542615 0.3652989 0.3542615 0.3541855 1.8710573 0.49706191 7.043097
    ##     4: 0.3195523 0.3188303 0.3195523 0.3186921 0.5412832 0.16162940 1.812712
    ##     5: 0.6080318 0.6084157 0.6080318 0.6077573 0.8151935 0.37337733 1.779809
    ##    ---                                                                      
    ## 22211: 0.3983666 0.3949865 0.3983666 0.3981244 0.6620237 0.25419512 1.724169
    ## 22212: 0.9016133 0.9015048 0.9016133 0.9016144 1.0493838 0.48869230 2.253373
    ## 22213: 0.8691759 0.8691994 0.8691759 0.8691733 1.0561974 0.55103934 2.024453
    ## 22214: 0.4496526 0.4478541 0.4496526 0.4498007 0.5245510 0.09847349 2.794191
    ## 22215: 0.9998490 0.9998490 0.9998490 0.9998490 0.9998604 0.23547784 4.245498
    ##         P.adjust LRT.adjust Wald.adjust LogRank.adjust
    ##     1: 0.9999969  0.9999969   0.9999969      0.9999969
    ##     2: 0.9999969  0.9999969   0.9999969      0.9999969
    ##     3: 0.9999969  0.9999969   0.9999969      0.9999969
    ##     4: 0.9999969  0.9999969   0.9999969      0.9999969
    ##     5: 0.9999969  0.9999969   0.9999969      0.9999969
    ##    ---                                                
    ## 22211: 0.9999969  0.9999969   0.9999969      0.9999969
    ## 22212: 0.9999969  0.9999969   0.9999969      0.9999969
    ## 22213: 0.9999969  0.9999969   0.9999969      0.9999969
    ## 22214: 0.9999969  0.9999969   0.9999969      0.9999969
    ## 22215: 0.9999969  0.9999969   0.9999969      0.9999969

We now take the top probes from the model by Log Rank p-value and use
<i>biomaRt</i> to look up the corresponding gene symbols.

*not run*

``` r
  res5 <- res5[order(res5$LogRank, decreasing = FALSE),]
  final <- subset(res5, LogRank < 0.01)
  probes <- gsub('^X', '', final$Variable)
  library(biomaRt)
  mart <- useMart('ENSEMBL_MART_ENSEMBL', host = 'useast.ensembl.org')
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  annotLookup <- getBM(mart = mart,
    attributes = c('affy_hg_u133a',
      'ensembl_gene_id',
      'gene_biotype',
      'external_gene_name'),
    filter = 'affy_hg_u133a',
    values = probes,
    uniqueRows = TRUE)
```

Two of the top hits include <i>CXCL12</i> and <i>MMP10</i>. High
expression of <i>CXCL12</i> was previously associated with good
progression free and overall survival in breast cancer in (doi:
10.1016/j.cca.2018.05.041.)\[<https://www.ncbi.nlm.nih.gov/pubmed/29800557>\]
, whilst high expression of <i>MMP10</i> was associated with poor
prognosis in colon cancer in (doi:
10.1186/s12885-016-2515-7)\[<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4950722/>\].

We can further explore the role of these genes to RFS by dividing their
gene expression Z-scores into tertiles for low, mid, and high
expression:

``` r
  # extract RFS and probe data for downstream analysis
  survplotdata <- coxdata[,c('Time.RFS', 'Distant.RFS',
    'X203666_at', 'X205680_at')]
  colnames(survplotdata) <- c('Time.RFS', 'Distant.RFS',
    'CXCL12', 'MMP10')
  # set Z-scale cut-offs for high and low expression
  highExpr <- 1.0
  lowExpr <- 1.0
  # encode the expression for CXCL12 and MMP10 as low, mid, and high
  survplotdata$CXCL12 <- ifelse(survplotdata$CXCL12 >= highExpr, 'High',
    ifelse(x <= lowExpr, 'Low', 'Mid'))
  survplotdata$MMP10 <- ifelse(survplotdata$MMP10 >= highExpr, 'High',
    ifelse(x <= lowExpr, 'Low', 'Mid'))
  # relevel the factors to have mid as the reference level
  survplotdata$CXCL12 <- factor(survplotdata$CXCL12,
    levels = c('Mid', 'Low', 'High'))
  survplotdata$MMP10 <- factor(survplotdata$MMP10,
    levels = c('Mid', 'Low', 'High'))
```

Plot the survival curves and place Log Rank p-value in the plots:

``` r
  library(survminer)
  ggsurvplot(survfit(Surv(Time.RFS, Distant.RFS) ~ CXCL12,
    data = survplotdata),
    data = survplotdata,
    risk.table = TRUE,
    pval = TRUE,
    break.time.by = 500,
    ggtheme = theme_minimal(),
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE)
```

![Survival analysis via Cox Proportional Hazards
regression.](README_files/figure-gfm/coxphParallel5-1.png)

``` r
  ggsurvplot(survfit(Surv(Time.RFS, Distant.RFS) ~ MMP10,
    data = survplotdata),
    data = survplotdata,
    risk.table = TRUE,
    pval = TRUE,
    break.time.by = 500,
    ggtheme = theme_minimal(),
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE)
```

![Survival analysis via Cox Proportional Hazards
regression.](README_files/figure-gfm/coxphParallel5-2.png)

## Perform a conditional logistic regression

In this example, we will re-use the Cox data for the purpose of
performing conditional logistic regression with tumour grade as our
grouping / matching factor. For this example, we will use ER status as
the dependent variable and also adjust for age.

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
```

    ##        Variable         Term       Beta StandardError         Z           P
    ## 1:   X204667_at   X204667_at  0.9940504     0.3628087  2.739875 0.006146252
    ## 2:   X205225_at   X205225_at  0.4444556     0.1633857  2.720285 0.006522559
    ## 3: X207813_s_at X207813_s_at  0.8218501     0.3050777  2.693904 0.007062046
    ## 4:   X212108_at   X212108_at  1.9610211     0.7607284  2.577820 0.009942574
    ## 5: X219497_s_at X219497_s_at -1.0249671     0.3541401 -2.894242 0.003800756
    ##            LRT       Wald    LogRank        HR   HRlower   HRupper
    ## 1: 0.006808415 0.02212540 0.02104525 2.7021573 1.3270501  5.502169
    ## 2: 0.010783544 0.01941078 0.01701248 1.5596409 1.1322713  2.148319
    ## 3: 0.037459927 0.02449358 0.02424809 2.2747043 1.2509569  4.136257
    ## 4: 0.033447973 0.03356050 0.03384960 7.1065797 1.6000274 31.564132
    ## 5: 0.005153233 0.01387183 0.01183245 0.3588083 0.1792329  0.718302

*not run*

``` r
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

Oestrogen receptor (<i>ESR1</i>) comes out - makes sense\! Also,
although 204667\_at is not listed in <i>biomaRt</i>, it overlaps an exon
of <i>FOXA1</i>, which also makes sense in relation to oestrogen
signalling.

    ##            used  (Mb) gc trigger  (Mb) max used  (Mb)
    ## Ncells  7061891 377.2   11504125 614.4 11504125 614.4
    ## Vcells 28488913 217.4   72389158 552.3 72389120 552.3

# Advanced features

Advanced features include the ability to modify block size, choose
different numbers of cores, enable ‘nested’ parallel processing, modify
limits for confidence intervals, and exclude certain model terms from
output.

## Speed up processing

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

With 2 cores instead of the default of 3, coupled with nestedParallel
being enabled, a total of 2 x 2 = 4 threads will be used.

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
    ##  11.352   4.656  10.412

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
    ##   0.988   0.164  10.084

Focusing on the elapsed time (as system time only reports time from the
last core that finished), we can see that nested processing has
negligible improvement or may actually be slower under certain
conditions when tested over a small number of variables. This is likely
due to the system being slowed by simply managing the larger number of
threads. Nested processing’s benefits can only be gained when processing
a large number of variables:

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
    ## 306.584  31.132 206.879

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
    ## 262.116   1.616 265.749

Performance is system-dependent and even increasing cores may not result
in huge gains in time. Performance is a trade-off between cores, forked
threads, blocksize, and the number of terms in each model.

### ~40000 tests; blocksize, 5000; cores, 3; nestedParallel, TRUE

In this example, we choose a large blocksize and 3 cores. With
nestedParallel enabled, this translates to 9 simultaneous threads.

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
    ## 616.236  26.196 309.488

## Modify confidence intervals

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

    ##       Variable         Term       Beta StandardError        Z        P       OR
    ##    1:    gene1        gene1 0.06333797     0.0546336 1.159323 0.246325  1.06539
    ##    2:    gene1 cellB:dosage 1.46710693     1.1706922 1.253196 0.210134  4.33667
    ##    3:    gene1 cellT:dosage 0.75012944     1.0646934 0.704550 0.481091  2.11727
    ##    4:    gene2        gene2 0.08600691     0.0587386 1.464231 0.143131  1.08981
    ##    5:    gene2 cellB:dosage 1.42579282     1.1307990 1.260872 0.207355  4.16116
    ##   ---                                                                          
    ## 1487:  gene496 cellB:dosage 2.70464635     1.7641808 1.533089 0.125254 14.94903
    ## 1488:  gene496 cellT:dosage 0.72487045     1.1111015 0.652389 0.514150  2.06446
    ## 1489:  gene497      gene497 0.00837784     0.0547597 0.152993 0.878404  1.00841
    ## 1490:  gene497 cellB:dosage 1.14481873     0.9726219 1.177044 0.239178  3.14187
    ## 1491:  gene497 cellT:dosage 0.48628896     0.9857062 0.493341 0.621772  1.62627
    ##        ORlower    ORupper
    ##    1: 0.925530    1.22638
    ##    2: 0.212589   88.46529
    ##    3: 0.136376   32.87124
    ##    4: 0.936792    1.26783
    ##    5: 0.226061   76.59547
    ##   ---                    
    ## 1487: 0.158884 1406.52161
    ## 1488: 0.117992   36.12114
    ## 1489: 0.875751    1.16117
    ## 1490: 0.256535   38.47954
    ## 1491: 0.128385   20.60018

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

    ##       Variable         Term       Beta StandardError        Z        P       OR
    ##    1:    gene1        gene1 0.06333797     0.0546336 1.159323 0.246325  1.06539
    ##    2:    gene1 cellB:dosage 1.46710693     1.1706922 1.253196 0.210134  4.33667
    ##    3:    gene1 cellT:dosage 0.75012944     1.0646934 0.704550 0.481091  2.11727
    ##    4:    gene2        gene2 0.08600691     0.0587386 1.464231 0.143131  1.08981
    ##    5:    gene2 cellB:dosage 1.42579282     1.1307990 1.260872 0.207355  4.16116
    ##   ---                                                                          
    ## 1487:  gene496 cellB:dosage 2.70464635     1.7641808 1.533089 0.125254 14.94903
    ## 1488:  gene496 cellT:dosage 0.72487045     1.1111015 0.652389 0.514150  2.06446
    ## 1489:  gene497      gene497 0.00837784     0.0547597 0.152993 0.878404  1.00841
    ## 1490:  gene497 cellB:dosage 1.14481873     0.9726219 1.177044 0.239178  3.14187
    ## 1491:  gene497 cellT:dosage 0.48628896     0.9857062 0.493341 0.621772  1.62627
    ##        ORlower   ORupper
    ##    1: 0.957201   1.18580
    ##    2: 0.437181  43.01812
    ##    3: 0.262729  17.06262
    ##    4: 0.971301   1.22279
    ##    5: 0.453603  38.17260
    ##   ---                   
    ## 1487: 0.470912 474.55485
    ## 1488: 0.233903  18.22127
    ## 1489: 0.905789   1.12266
    ## 1490: 0.466972  21.13906
    ## 1491: 0.235591  11.22606

## Remove some terms from output / include the intercept

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

    ##      Variable        Term        Beta StandardError         Z        P       OR
    ##   1:    gene1 (Intercept) -1.35754580     1.0199412 -1.331004 0.183188 0.257291
    ##   2:    gene1       gene1  0.06333797     0.0546336  1.159323 0.246325 1.065387
    ##   3:    gene2 (Intercept) -1.92893568     1.2193210 -1.581975 0.113655 0.145303
    ##   4:    gene2       gene2  0.08600691     0.0587386  1.464231 0.143131 1.089814
    ##   5:    gene3 (Intercept) -1.36553647     0.9199838 -1.484305 0.137728 0.255244
    ##  ---                                                                           
    ## 990:  gene495     gene495  0.02963303     0.0426166  0.695340 0.486842 1.030076
    ## 991:  gene496 (Intercept) -0.02762196     0.8035955 -0.034373 0.972580 0.972756
    ## 992:  gene496     gene496 -0.10668807     0.0699128 -1.526017 0.127006 0.898806
    ## 993:  gene497 (Intercept) -0.66652515     0.8262057 -0.806730 0.419822 0.513490
    ## 994:  gene497     gene497  0.00837784     0.0547597  0.152993 0.878404 1.008413
    ##        ORlower ORupper
    ##   1: 0.0348538 1.89933
    ##   2: 0.9572010 1.18580
    ##   3: 0.0133164 1.58548
    ##   4: 0.9713012 1.22279
    ##   5: 0.0420594 1.54898
    ##  ---                  
    ## 990: 0.9475326 1.11981
    ## 991: 0.2013642 4.69922
    ## 992: 0.7837113 1.03080
    ## 993: 0.1016867 2.59298
    ## 994: 0.9057888 1.12266

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

    ##      Variable    Term        Beta StandardError         Z        P       OR
    ##   1:    gene1   gene1  0.06333797     0.0546336  1.159323 0.246325 1.065387
    ##   2:    gene2   gene2  0.08600691     0.0587386  1.464231 0.143131 1.089814
    ##   3:    gene3   gene3  0.05882164     0.0412150  1.427189 0.153525 1.060586
    ##   4:    gene4   gene4  0.17931464     0.1044105  1.717401 0.085906 1.196397
    ##   5:    gene5   gene5  0.09056917     0.0639784  1.415622 0.156886 1.094797
    ##  ---                                                                       
    ## 493:  gene493 gene493 -0.06615597     0.0798339 -0.828670 0.407291 0.935985
    ## 494:  gene494 gene494  0.06632122     0.0630053  1.052630 0.292511 1.068570
    ## 495:  gene495 gene495  0.02963303     0.0426166  0.695340 0.486842 1.030076
    ## 496:  gene496 gene496 -0.10668807     0.0699128 -1.526017 0.127006 0.898806
    ## 497:  gene497 gene497  0.00837784     0.0547597  0.152993 0.878404 1.008413
    ##       ORlower ORupper
    ##   1: 0.957201 1.18580
    ##   2: 0.971301 1.22279
    ##   3: 0.978281 1.14982
    ##   4: 0.974992 1.46808
    ##   5: 0.965773 1.24106
    ##  ---                 
    ## 493: 0.800413 1.09452
    ## 494: 0.944436 1.20902
    ## 495: 0.947533 1.11981
    ## 496: 0.783711 1.03080
    ## 497: 0.905789 1.12266

# Acknowledgments

Thanks to Horácio Montenegro and GenoMax for testing cross-platform
differences, and Wolfgang Huber for providing the nudge that FDR
correction needed to be implemented.

Thanks to Michael Barnes in London for introducing me to parallel
processing in R.

Finally, thanks to Juan Celedón at Children’s Hospital of Pittsburgh.

Sarega Gurudas, whose suggestion led to the implementation of survey
weights via svyglm.

# Session info

``` r
sessionInfo()
```

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.7 LTS
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
    ##  [1] survminer_0.4.9             ggpubr_0.4.0               
    ##  [3] ggplot2_3.3.3               GEOquery_2.56.0            
    ##  [5] DESeq2_1.28.1               magrittr_2.0.1             
    ##  [7] airway_1.8.0                SummarizedExperiment_1.18.2
    ##  [9] DelayedArray_0.14.1         matrixStats_0.57.0         
    ## [11] Biobase_2.48.0              GenomicRanges_1.40.0       
    ## [13] GenomeInfoDb_1.24.2         IRanges_2.22.2             
    ## [15] S4Vectors_0.26.1            BiocGenerics_0.34.0        
    ## [17] RegParallel_1.11.1          arm_1.11-2                 
    ## [19] lme4_1.1-26                 Matrix_1.3-2               
    ## [21] MASS_7.3-53                 survival_3.2-7             
    ## [23] stringr_1.4.0               data.table_1.13.6          
    ## [25] doParallel_1.0.16           iterators_1.0.13           
    ## [27] foreach_1.5.1               knitr_1.31                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1           backports_1.2.1        Hmisc_4.4-2           
    ##   [4] splines_4.0.3          BiocParallel_1.22.0    digest_0.6.27         
    ##   [7] htmltools_0.5.1.1      fansi_0.4.2            checkmate_2.0.0       
    ##  [10] memoise_2.0.0          cluster_2.1.0          openxlsx_4.2.3        
    ##  [13] limma_3.44.3           readr_1.4.0            annotate_1.66.0       
    ##  [16] jpeg_0.1-8.1           colorspace_2.0-0       blob_1.2.1            
    ##  [19] haven_2.3.1            xfun_0.20              dplyr_1.0.3           
    ##  [22] crayon_1.3.4           RCurl_1.98-1.2         genefilter_1.70.0     
    ##  [25] zoo_1.8-8              glue_1.4.2             gtable_0.3.0          
    ##  [28] zlibbioc_1.34.0        XVector_0.28.0         car_3.0-10            
    ##  [31] abind_1.4-5            scales_1.1.1           DBI_1.1.1             
    ##  [34] rstatix_0.7.0          Rcpp_1.0.6             xtable_1.8-4          
    ##  [37] gridtext_0.1.4         htmlTable_2.1.0        foreign_0.8-81        
    ##  [40] bit_4.0.4              km.ci_0.5-2            Formula_1.2-4         
    ##  [43] htmlwidgets_1.5.3      RColorBrewer_1.1-2     ellipsis_0.3.1        
    ##  [46] pkgconfig_2.0.3        XML_3.99-0.5           farver_2.0.3          
    ##  [49] nnet_7.3-15            locfit_1.5-9.4         tidyselect_1.1.0      
    ##  [52] labeling_0.4.2         rlang_0.4.10           AnnotationDbi_1.53.0  
    ##  [55] munsell_0.5.0          cellranger_1.1.0       tools_4.0.3           
    ##  [58] cachem_1.0.1           cli_2.2.0              generics_0.1.0        
    ##  [61] RSQLite_2.2.3          broom_0.7.5            evaluate_0.14         
    ##  [64] fastmap_1.1.0          yaml_2.2.1             bit64_4.0.5           
    ##  [67] zip_2.1.1              survMisc_0.5.5         purrr_0.3.4           
    ##  [70] nlme_3.1-151           xml2_1.3.2             compiler_4.0.3        
    ##  [73] rstudioapi_0.13        curl_4.3               png_0.1-7             
    ##  [76] ggsignif_0.6.0         tibble_3.0.1           statmod_1.4.35        
    ##  [79] geneplotter_1.66.0     stringi_1.5.3          highr_0.8             
    ##  [82] forcats_0.5.1          lattice_0.20-41        markdown_1.1          
    ##  [85] nloptr_1.2.2.2         KMsurv_0.1-5           vctrs_0.3.6           
    ##  [88] pillar_1.4.7           lifecycle_0.2.0        bitops_1.0-6          
    ##  [91] R6_2.5.0               latticeExtra_0.6-29    gridExtra_2.3         
    ##  [94] rio_0.5.16             codetools_0.2-18       boot_1.3-26           
    ##  [97] assertthat_0.2.1       withr_2.4.1            GenomeInfoDbData_1.2.3
    ## [100] hms_1.0.0              ggtext_0.1.1           grid_4.0.3            
    ## [103] rpart_4.1-15           tidyr_1.1.2            coda_0.19-4           
    ## [106] minqa_1.2.4            rmarkdown_2.6          carData_3.0-4         
    ## [109] base64enc_0.1-3

# References

Blighe and Lasky-Su (2018)

<div id="refs" class="references">

<div id="ref-RegParallel">

Blighe, K, and J Lasky-Su. 2018. “RegParallel: Standard regression
functions in R enabled for parallel processing over large data-frames.”
<https://github.com/kevinblighe/RegParallel.>

</div>

</div>
