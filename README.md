Standard regression functions in R enabled for parallel processing over large data-frames
================
Kevin Blighe
2018-09-25

-   [Introduction.](#introduction.)
-   [Installation.](#installation.)
    -   [1. Download the package from Bioconductor.](#download-the-package-from-bioconductor.)
    -   [2. Load the package into R session.](#load-the-package-into-r-session.)
-   [Quick start.](#quick-start.)
    -   [Perform the most basic regression analysis.](#perform-the-most-basic-regression-analysis.)
    -   [Perform a basic linear regression.](#perform-a-basic-linear-regression.)
    -   [Perform the most basic negative binomial logistic regression analysis with glm.nbParallel:](#perform-the-most-basic-negative-binomial-logistic-regression-analysis-with-glm.nbparallel)
    -   [Survival analysis via Cox Proportional Hazards regression.](#survival-analysis-via-cox-proportional-hazards-regression.)
    -   [Perform a conditional logistic regression](#perform-a-conditional-logistic-regression)
-   [Advanced features.](#advanced-features.)
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

Introduction.
=============

In many analyses, a large amount of variables have to be tested independently against the trait/endpoint of interest, and also adjusted for covariates and confounding factors at the same time. The major botteleneck in these is the amount of time that it takes to complete these analyses.

With <i>RegParallel</i>, any number of tests can be performed simultaneously. On a 12-core system, 144 variables can be tested simultaneously, with 1000s of variables processed in a matter of seconds via 'nested' parallel processing.

Works for logistic regression, linear regression, conditional logistic regression, Cox proportional hazards and survival models, Bayesian logistic regression, and negative binomial regression.

Installation.
=============

1. Download the package from Bioconductor.
------------------------------------------

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) 
install.packages("BiocManager")

BiocManager::install("RegParallel")
```

Note: to install development version:

``` r
devtools::install_github("kevinblighe/RegParallel")
```

2. Load the package into R session.
-----------------------------------

``` r
library(RegParallel)
```

Quick start.
============

For this quick start, we will follow the tutorial (from Section 3.1) of [RNA-seq workflow: gene-level exploratory analysis and differential expression](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html). Specifically, we will load the 'airway' data, where different airway smooth muscle cells were treated with dexamethasone.

``` r
library(airway)

library(magrittr)

data("airway")

airway$dex %<>% relevel("untrt")
```

Normalise the raw counts in DESeq2 and produce regularised log counts:

``` r
library(DESeq2)

dds <- DESeqDataSet(airway, design = ~dex + cell)
dds <- DESeq(dds, betaPrior = FALSE)
rlogcounts <- assay(rlog(dds, blind = FALSE))
rlogdata <- data.frame(colData(airway), t(rlogcounts))
```

Perform the most basic regression analysis.
-------------------------------------------

Here, we fit a binomial logistic regression model to the data via glmParallel:

``` r
res1 <- RegParallel(data = rlogdata[, 1:3000], formula = "dex ~ [*]", FUN = function(formula, 
    data) glm(formula = formula, data = data, family = binomial(link = "logit")), 
    FUNtype = "glm", variables = colnames(rlogdata)[10:3000])

res1[order(res1$P, decreasing = FALSE), ]
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
    ## 2821: -0.0004028919 0.99967854  3.514359e-95 0.000000e+00           NA

Perform a basic linear regression.
----------------------------------

Here, we will perform the linear regression using both glmParallel and lmParallel. We will appreciate that a linear regression is the same using either function with the default settings.

Regularised log counts from our DESeq2 data will be used.

``` r
rlogdata <- rlogdata[, 1:2000]

res2 <- RegParallel(data = rlogdata, formula = "[*] ~ cell", FUN = function(formula, 
    data) glm(formula = formula, data = data, method = "glm.fit"), FUNtype = "glm", 
    variables = colnames(rlogdata)[10:ncol(rlogdata)])

res3 <- RegParallel(data = rlogdata, formula = "[*] ~ cell", FUN = function(formula, 
    data) lm(formula = formula, data = data), FUNtype = "lm", variables = colnames(rlogdata)[10:ncol(rlogdata)])

subset(res2, P < 0.05)
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
subset(res3, P < 0.05)
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

Perform the most basic negative binomial logistic regression analysis with glm.nbParallel:
------------------------------------------------------------------------------------------

Here, we will utilise normalised, unlogged counts from DESeq2.

``` r
nbcounts <- round(counts(dds, normalized = TRUE), 0)
nbdata <- data.frame(colData(airway), t(nbcounts))

res4 <- RegParallel(data = nbdata[, 1:3000], formula = "[*] ~ dex", FUN = function(formula, 
    data) glm.nb(formula = formula, data = data), FUNtype = "glm.nb", variables = colnames(nbdata)[10:3000])

res4[order(res4$Theta, decreasing = TRUE), ]
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
    ##              OR    ORlower     ORupper
    ##    1: 0.8699790 0.83880014   0.9023169
    ##    2: 0.9701750 0.94055425   1.0007286
    ##    3: 0.9971811 0.94658900   1.0504772
    ##    4: 0.8796624 0.83904459   0.9222465
    ##    5: 1.1437649 1.08420938   1.2065918
    ##   ---                                 
    ## 2817: 0.2500000 0.00834686   7.4878457
    ## 2818: 4.0000000 0.13354976 119.8055313
    ## 2819: 0.2500000 0.00834686   7.4878457
    ## 2820: 0.7272727 0.04715465  11.2168280
    ## 2821: 1.7368421 0.07355573  41.0113595

Survival analysis via Cox Proportional Hazards regression.
----------------------------------------------------------

For this example, we will load breast cancer gene expression data with recurrence free survival (RFS) from [Gene Expression Profiling in Breast Cancer: Understanding the Molecular Basis of Histologic Grade To Improve Prognosis](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2990). Specifically, we will encode each gene's expression into Low|Mid|High based on Z-scores and compare these against RFS while adjusting for tumour grade in a Cox Proportional Hazards model.

First, let's read in and prepare the data

``` r
library(Biobase)

library(GEOquery)

# load series and platform data from GEO
gset <- getGEO("GSE2990", GSEMatrix = TRUE, getGPL = FALSE)

x <- exprs(gset[[1]])

# remove Affymetrix control probes
x <- x[-grep("^AFFX", rownames(x)), ]

# transform the expression data to Z scores
x <- t(scale(t(x)))

# extract information of interest from the phenotype data (pdata)
idx <- which(colnames(pData(gset[[1]])) %in% c("age:ch1", "distant rfs:ch1", 
    "er:ch1", "ggi:ch1", "grade:ch1", "node:ch1", "size:ch1", "time rfs:ch1"))

metadata <- data.frame(pData(gset[[1]])[, idx], row.names = rownames(pData(gset[[1]])))

# remove samples from the pdata that have any NA value
discard <- apply(metadata, 1, function(x) any(is.na(x)))
metadata <- metadata[!discard, ]

# filter the Z-scores expression data to match the samples in our pdata
x <- x[, which(colnames(x) %in% rownames(metadata))]

# check that sample names match exactly between pdata and Z-scores
all((colnames(x) == rownames(metadata)) == TRUE)
```

    ## [1] TRUE

``` r
# create a merged pdata and Z-scores object
coxdata <- data.frame(metadata, t(x))

# tidy column names
colnames(coxdata)[1:8] <- c("Age", "Distant.RFS", "ER", "GGI", "Grade", "Node", 
    "Size", "Time.RFS")

# prepare certain phenotypes
coxdata$Age <- as.numeric(gsub("^KJ", "", coxdata$Age))
coxdata$Distant.RFS <- as.numeric(coxdata$Distant.RFS)
coxdata$ER <- factor(coxdata$ER, levels = c(0, 1))
coxdata$Grade <- factor(coxdata$Grade, levels = c(1, 2, 3))
coxdata$Time.RFS <- as.numeric(gsub("^KJX|^KJ", "", coxdata$Time.RFS))
```

With the data prepared, we can now apply a Cox Proportional Hazards model indpendently for each probe in the dataset against RFS.

In this we also increase the default blocksize to 2000 in order to speed up the analysis.

``` r
library(survival)

res5 <- RegParallel(data = coxdata, formula = "Surv(Time.RFS, Distant.RFS) ~ [*]", 
    FUN = function(formula, data) coxph(formula = formula, data = data, ties = "breslow", 
        singular.ok = TRUE), FUNtype = "coxph", variables = colnames(coxdata)[9:ncol(coxdata)], 
    blocksize = 2000)

res5 <- res5[!is.na(res5$P), ]

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
    ##                P       LRT      Wald   LogRank        HR    HRlower
    ##     1: 0.2848529 0.2826716 0.2848529 0.2848400 1.4594563 0.72994385
    ##     2: 0.6047873 0.6085603 0.6047873 0.6046839 1.1249515 0.72024775
    ##     3: 0.3542615 0.3652989 0.3542615 0.3541855 1.8710573 0.49706191
    ##     4: 0.3195523 0.3188303 0.3195523 0.3186921 0.5412832 0.16162940
    ##     5: 0.6080318 0.6084157 0.6080318 0.6077573 0.8151935 0.37337733
    ##    ---                                                             
    ## 22211: 0.3983666 0.3949865 0.3983666 0.3981244 0.6620237 0.25419512
    ## 22212: 0.9016133 0.9015048 0.9016133 0.9016144 1.0493838 0.48869230
    ## 22213: 0.8691759 0.8691994 0.8691759 0.8691733 1.0561974 0.55103934
    ## 22214: 0.4496526 0.4478541 0.4496526 0.4498007 0.5245510 0.09847349
    ## 22215: 0.9998490 0.9998490 0.9998490 0.9998490 0.9998604 0.23547784
    ##         HRupper
    ##     1: 2.918050
    ##     2: 1.757056
    ##     3: 7.043097
    ##     4: 1.812712
    ##     5: 1.779809
    ##    ---         
    ## 22211: 1.724169
    ## 22212: 2.253373
    ## 22213: 2.024453
    ## 22214: 2.794191
    ## 22215: 4.245498

We now take the top probes from the model by Log Rank p-value and use <i>biomaRt</i> to look up the corresponding gene symbols.

``` r
res5 <- res5[order(res5$LogRank, decreasing = FALSE), ]

final <- subset(res5, LogRank < 0.01)

probes <- gsub("^X", "", final$Variable)

library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(mart = mart, attributes = c("affy_hg_u133a", "ensembl_gene_id", 
    "gene_biotype", "external_gene_name"), filter = "affy_hg_u133a", values = probes, 
    uniqueRows = TRUE)
```

Two of the top hits include <i>CXCL12</i> and <i>MMP10</i>. High expression of <i>CXCL12</i> was previously associated with good progression free and overall survival in breast cancer in (doi: 10.1016/j.cca.2018.05.041.)\[<https://www.ncbi.nlm.nih.gov/pubmed/29800557>\], whilst high expression of <i>MMP10</i> was associated with poor prognosis in colon cancer in (doi: 10.1186/s12885-016-2515-7)\[<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4950722/>\].

We can further explore the role of these genes to RFS by dividing their gene expression Z-scores into tertiles for low, mid, and high expression:

``` r
# extract RFS and probe data for downstream analysis
survplotdata <- coxdata[, c("Time.RFS", "Distant.RFS", "X203666_at", "X205680_at")]
colnames(survplotdata) <- c("Time.RFS", "Distant.RFS", "CXCL12", "MMP10")

# set Z-scale cut-offs for high and low expression
highExpr <- 1
lowExpr <- 1

# encode the expression for CXCL12 and MMP10 as low, mid, and high
survplotdata$CXCL12 <- ifelse(survplotdata$CXCL12 >= highExpr, "High", ifelse(x <= 
    lowExpr, "Low", "Mid"))
survplotdata$MMP10 <- ifelse(survplotdata$MMP10 >= highExpr, "High", ifelse(x <= 
    lowExpr, "Low", "Mid"))

# relevel the factors to have mid as the reference level
survplotdata$CXCL12 <- factor(survplotdata$CXCL12, levels = c("Mid", "Low", 
    "High"))
survplotdata$MMP10 <- factor(survplotdata$MMP10, levels = c("Mid", "Low", "High"))
```

Plot the survival curves and place Log RAnk p-value in the plots:

``` r
library(survminer)

ggsurvplot(survfit(Surv(Time.RFS, Distant.RFS) ~ CXCL12, data = survplotdata), 
    data = survplotdata, risk.table = TRUE, pval = TRUE, break.time.by = 500, 
    ggtheme = theme_minimal(), risk.table.y.text.col = TRUE, risk.table.y.text = FALSE)
```

![Survival analysis via Cox Proportional Hazards regression.](README_files/figure-markdown_github/coxphParallel5-1.png)

``` r
ggsurvplot(survfit(Surv(Time.RFS, Distant.RFS) ~ MMP10, data = survplotdata), 
    data = survplotdata, risk.table = TRUE, pval = TRUE, break.time.by = 500, 
    ggtheme = theme_minimal(), risk.table.y.text.col = TRUE, risk.table.y.text = FALSE)
```

![Survival analysis via Cox Proportional Hazards regression.](README_files/figure-markdown_github/coxphParallel5-2.png)

Perform a conditional logistic regression
-----------------------------------------

In this example, we will re-use the Cox data for the purpose of performing conditional logistic regression with tumour grade as our grouping / matching factor. For this example, we will use ER status as the dependent variable and also adjust for age.

``` r
x <- exprs(gset[[1]])
x <- x[-grep("^AFFX", rownames(x)), ]
x <- scale(x)
x <- x[, which(colnames(x) %in% rownames(metadata))]
coxdata <- data.frame(metadata, t(x))
colnames(coxdata)[1:8] <- c("Age", "Distant.RFS", "ER", "GGI", "Grade", "Node", 
    "Size", "Time.RFS")
coxdata$Age <- as.numeric(gsub("^KJ", "", coxdata$Age))
coxdata$Grade <- factor(coxdata$Grade, levels = c(1, 2, 3))
coxdata$ER <- as.numeric(coxdata$ER)
coxdata <- coxdata[!is.na(coxdata$ER), ]

res6 <- RegParallel(data = coxdata, formula = "ER ~ [*] + Age + strata(Grade)", 
    FUN = function(formula, data) clogit(formula = formula, data = data, method = "breslow"), 
    FUNtype = "clogit", variables = colnames(coxdata)[9:ncol(coxdata)], blocksize = 2000)

subset(res6, P < 0.01)
```

    ##        Variable         Term       Beta StandardError         Z
    ## 1:   X204667_at   X204667_at  0.9940504     0.3628087  2.739875
    ## 2:   X205225_at   X205225_at  0.4444556     0.1633857  2.720285
    ## 3: X207813_s_at X207813_s_at  0.8218501     0.3050777  2.693904
    ## 4:   X212108_at   X212108_at  1.9610211     0.7607284  2.577820
    ## 5: X219497_s_at X219497_s_at -1.0249671     0.3541401 -2.894242
    ##              P         LRT       Wald    LogRank        HR   HRlower
    ## 1: 0.006146252 0.006808415 0.02212540 0.02104525 2.7021573 1.3270501
    ## 2: 0.006522559 0.010783544 0.01941078 0.01701248 1.5596409 1.1322713
    ## 3: 0.007062046 0.037459927 0.02449358 0.02424809 2.2747043 1.2509569
    ## 4: 0.009942574 0.033447973 0.03356050 0.03384960 7.1065797 1.6000274
    ## 5: 0.003800756 0.005153233 0.01387183 0.01183245 0.3588083 0.1792329
    ##      HRupper
    ## 1:  5.502169
    ## 2:  2.148319
    ## 3:  4.136257
    ## 4: 31.564132
    ## 5:  0.718302

``` r
getBM(mart = mart, attributes = c("affy_hg_u133a", "ensembl_gene_id", "gene_biotype", 
    "external_gene_name"), filter = "affy_hg_u133a", values = c("204667_at", 
    "205225_at", "207813_s_at", "212108_at", "219497_s_at"), uniqueRows = TRUE)
```

    ##   affy_hg_u133a ensembl_gene_id   gene_biotype external_gene_name
    ## 1     205225_at ENSG00000091831 protein_coding               ESR1
    ## 2     212108_at ENSG00000113194 protein_coding               FAF2
    ## 3   207813_s_at ENSG00000161513 protein_coding               FDXR
    ## 4   219497_s_at ENSG00000119866 protein_coding             BCL11A

<i>Oestrogen receptor (ESR1)</i> comes out top - makes sense! Also, although 204667\_at is not listed in <i>biomaRt</i>, it overlaps an exon of <i>FOXA1</i>, which also makes sense in relation to oestrogen signalling.

Advanced features.
==================

Advanced features include the ability to modify block size, choose different numbers of cores, enable 'nested' parallel processing, modify limits for confidence intervals, and exclude certain model terms from output.

Speed up processing
-------------------

First create some test data for the purpose of benchmarking:

``` r
options(scipen = 10)
options(digits = 6)

col <- 60000
row <- 20
mat <- matrix(rexp(col * row, rate = 0.1), ncol = col)
colnames(mat) <- paste0("gene", 1:ncol(mat))
rownames(mat) <- paste0("sample", 1:nrow(mat))

modelling <- data.frame(cell = rep(c("B", "T"), nrow(mat)/2), group = c(rep(c("treatment"), 
    nrow(mat)/2), rep(c("control"), nrow(mat)/2)), dosage = t(data.frame(matrix(rexp(row, 
    rate = 1), ncol = row))), mat, row.names = rownames(mat))
```

### ~2000 tests; blocksize, 500; cores, 2; nestedParallel, TRUE

``` r
df <- modelling[, 1:2000]
variables <- colnames(df)[4:ncol(df)]

ptm <- proc.time()
res <- RegParallel(data = df, formula = "factor(group) ~ [*] + (cell:dosage) ^ 2", 
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = "logit"), 
        method = "glm.fit"), FUNtype = "glm", variables = variables, blocksize = 500, 
    cores = 2, nestedParallel = TRUE)
proc.time() - ptm
```

    ##    user  system elapsed 
    ##  12.876   4.764  11.599

### ~2000 tests; blocksize, 500; cores, 2; nestedParallel, FALSE

``` r
df <- modelling[, 1:2000]
variables <- colnames(df)[4:ncol(df)]

ptm <- proc.time()
res <- RegParallel(data = df, formula = "factor(group) ~ [*] + (cell:dosage) ^ 2", 
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = "logit"), 
        method = "glm.fit"), FUNtype = "glm", variables = variables, blocksize = 500, 
    cores = 2, nestedParallel = FALSE)
proc.time() - ptm
```

    ##    user  system elapsed 
    ##  10.660   2.340  11.637

Focusing on the elapsed time (as system time only reports time from the last core that finished), we can see that nested processing has no effect or is actually slower over a small number of variables. This is likely due to the system being slowed by simply managing the larger number of threads. Nested processing's benefits can only be gained when processing a large number of variables:

### ~40000 tests; blocksize, 2000; cores, 2; nestedParallel, TRUE

``` r
df <- modelling[, 1:40000]
variables <- colnames(df)[4:ncol(df)]

system.time(RegParallel(data = df, formula = "factor(group) ~ [*] + (cell:dosage) ^ 2", 
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = "logit"), 
        method = "glm.fit"), FUNtype = "glm", variables = variables, blocksize = 2000, 
    cores = 2, nestedParallel = TRUE))
```

    ##    user  system elapsed 
    ## 309.852  29.056 201.594

### ~40000 tests; blocksize, 2000; cores, 2; nestedParallel, FALSE

``` r
df <- modelling[, 1:40000]
variables <- colnames(df)[4:ncol(df)]

system.time(RegParallel(data = df, formula = "factor(group) ~ [*] + (cell:dosage) ^ 2", 
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = "logit"), 
        method = "glm.fit"), FUNtype = "glm", variables = variables, blocksize = 2000, 
    cores = 2, nestedParallel = FALSE))
```

    ##    user  system elapsed 
    ## 270.164  19.040 263.979

Performance is system-dependent and even increasing cores may not result in huge gains in time. Performance is a trade-off between cores, forked threads, blocksize, and the number of terms in each model.

### ~40000 tests; blocksize, 5000; cores, 3; nestedParallel, TRUE

``` r
df <- modelling[, 1:40000]
variables <- colnames(df)[4:ncol(df)]

system.time(RegParallel(data = df, formula = "factor(group) ~ [*] + (cell:dosage) ^ 2", 
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = "logit"), 
        method = "glm.fit"), FUNtype = "glm", variables = variables, blocksize = 500, 
    cores = 3, nestedParallel = TRUE))
```

    ##    user  system elapsed 
    ## 474.700 211.424 208.880

Modify confidence intervals
---------------------------

``` r
df <- modelling[, 1:1000]
variables <- colnames(df)[4:ncol(df)]

RegParallel(data = df, formula = "factor(group) ~ [*] + (cell:dosage) ^ 2", 
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = "logit"), 
        method = "glm.fit"), FUNtype = "glm", variables = variables, blocksize = 500, 
    cores = 3, nestedParallel = TRUE, conflevel = 99)
```

    ##       Variable         Term        Beta StandardError          Z        P
    ##    1:    gene1        gene1  0.16554942     0.1049014  1.5781428 0.114533
    ##    2:    gene1 cellB:dosage -0.41148815     0.5609013 -0.7336196 0.463181
    ##    3:    gene1 cellT:dosage  0.64680114     0.9348007  0.6919134 0.488992
    ##    4:    gene2        gene2  0.00108501     0.0584285  0.0185699 0.985184
    ##    5:    gene2 cellB:dosage -0.29819946     0.4797355 -0.6215914 0.534211
    ##   ---                                                                    
    ## 2987:  gene996 cellB:dosage -0.29819786     0.4794249 -0.6219908 0.533948
    ## 2988:  gene996 cellT:dosage  0.24119659     0.7776380  0.3101656 0.756435
    ## 2989:  gene997      gene997 -0.03448364     0.1154558 -0.2986740 0.765189
    ## 2990:  gene997 cellB:dosage -0.26387140     0.4923161 -0.5359797 0.591973
    ## 2991:  gene997 cellT:dosage  0.38915157     0.9256075  0.4204283 0.674173
    ##             OR  ORlower  ORupper
    ##    1: 1.180041 0.900632  1.54613
    ##    2: 0.662663 0.156254  2.81032
    ##    3: 1.909423 0.171858 21.21457
    ##    4: 1.001086 0.861210  1.16368
    ##    5: 0.742153 0.215689  2.55363
    ##   ---                           
    ## 2987: 0.742154 0.215862  2.55160
    ## 2988: 1.272771 0.171724  9.43343
    ## 2989: 0.966104 0.717575  1.30071
    ## 2990: 0.768072 0.216104  2.72986
    ## 2991: 1.475728 0.136006 16.01232

``` r
RegParallel(data = df, formula = "factor(group) ~ [*] + (cell:dosage) ^ 2", 
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = "logit"), 
        method = "glm.fit"), FUNtype = "glm", variables = variables, blocksize = 500, 
    cores = 3, nestedParallel = TRUE, conflevel = 95)
```

    ##       Variable         Term        Beta StandardError          Z        P
    ##    1:    gene1        gene1  0.16554942     0.1049014  1.5781428 0.114533
    ##    2:    gene1 cellB:dosage -0.41148815     0.5609013 -0.7336196 0.463181
    ##    3:    gene1 cellT:dosage  0.64680114     0.9348007  0.6919134 0.488992
    ##    4:    gene2        gene2  0.00108501     0.0584285  0.0185699 0.985184
    ##    5:    gene2 cellB:dosage -0.29819946     0.4797355 -0.6215914 0.534211
    ##   ---                                                                    
    ## 2987:  gene996 cellB:dosage -0.29819786     0.4794249 -0.6219908 0.533948
    ## 2988:  gene996 cellT:dosage  0.24119659     0.7776380  0.3101656 0.756435
    ## 2989:  gene997      gene997 -0.03448364     0.1154558 -0.2986740 0.765189
    ## 2990:  gene997 cellB:dosage -0.26387140     0.4923161 -0.5359797 0.591973
    ## 2991:  gene997 cellT:dosage  0.38915157     0.9256075  0.4204283 0.674173
    ##             OR  ORlower  ORupper
    ##    1: 1.180041 0.960738  1.44940
    ##    2: 0.662663 0.220726  1.98945
    ##    3: 1.909423 0.305632 11.92905
    ##    4: 1.001086 0.892764  1.12255
    ##    5: 0.742153 0.289829  1.90040
    ##   ---                           
    ## 2987: 0.742154 0.290006  1.89925
    ## 2988: 1.272771 0.277219  5.84355
    ## 2989: 0.966104 0.770456  1.21143
    ## 2990: 0.768072 0.292646  2.01587
    ## 2991: 1.475728 0.240507  9.05493

Remove some terms from output / include the intercept
-----------------------------------------------------

``` r
RegParallel(data = df, formula = "factor(group) ~ [*] + (cell:dosage) ^ 2", 
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = "logit"), 
        method = "glm.fit"), FUNtype = "glm", variables = variables, blocksize = 500, 
    cores = 3, nestedParallel = TRUE, conflevel = 95, excludeTerms = c("cell", 
        "dosage"), excludeIntercept = FALSE)
```

    ##       Variable        Term        Beta StandardError          Z        P
    ##    1:    gene1 (Intercept) -1.68014257     1.3018498 -1.2905809 0.196849
    ##    2:    gene1       gene1  0.16554942     0.1049014  1.5781428 0.114533
    ##    3:    gene2 (Intercept)  0.08587741     0.9571946  0.0897178 0.928511
    ##    4:    gene2       gene2  0.00108501     0.0584285  0.0185699 0.985184
    ##    5:    gene3 (Intercept) -0.27822010     0.8720416 -0.3190445 0.749693
    ##   ---                                                                   
    ## 1990:  gene995     gene995 -0.05314197     0.0795063 -0.6683996 0.503879
    ## 1991:  gene996 (Intercept)  0.08673635     0.8930973  0.0971186 0.922632
    ## 1992:  gene996     gene996  0.00102571     0.0493882  0.0207682 0.983431
    ## 1993:  gene997 (Intercept)  0.13794381     0.7477693  0.1844737 0.853642
    ## 1994:  gene997     gene997 -0.03448364     0.1154558 -0.2986740 0.765189
    ##             OR   ORlower ORupper
    ##    1: 0.186347 0.0145274 2.39034
    ##    2: 1.180041 0.9607380 1.44940
    ##    3: 1.089673 0.1669285 7.11314
    ##    4: 1.001086 0.8927642 1.12255
    ##    5: 0.757130 0.1370526 4.18267
    ##   ---                           
    ## 1990: 0.948245 0.8114184 1.10815
    ## 1991: 1.090609 0.1894364 6.27877
    ## 1992: 1.001026 0.9086700 1.10277
    ## 1993: 1.147911 0.2650976 4.97062
    ## 1994: 0.966104 0.7704560 1.21143

``` r
RegParallel(data = df, formula = "factor(group) ~ [*] + (cell:dosage) ^ 2", 
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = "logit"), 
        method = "glm.fit"), FUNtype = "glm", variables = variables, blocksize = 500, 
    cores = 3, nestedParallel = TRUE, conflevel = 95, excludeTerms = c("cell", 
        "dosage"), excludeIntercept = TRUE)
```

    ##      Variable    Term         Beta StandardError          Z        P
    ##   1:    gene1   gene1  0.165549419     0.1049014  1.5781428 0.114533
    ##   2:    gene2   gene2  0.001085014     0.0584285  0.0185699 0.985184
    ##   3:    gene3   gene3  0.075806046     0.0751546  1.0086683 0.313134
    ##   4:    gene4   gene4  0.061573792     0.0634040  0.9711345 0.331481
    ##   5:    gene5   gene5  0.000378096     0.0453364  0.0083398 0.993346
    ##  ---                                                                
    ## 993:  gene993 gene993 -0.026328218     0.0384402 -0.6849135 0.493399
    ## 994:  gene994 gene994  0.024994814     0.0475125  0.5260678 0.598841
    ## 995:  gene995 gene995 -0.053141965     0.0795063 -0.6683996 0.503879
    ## 996:  gene996 gene996  0.001025705     0.0493882  0.0207682 0.983431
    ## 997:  gene997 gene997 -0.034483637     0.1154558 -0.2986740 0.765189
    ##            OR  ORlower ORupper
    ##   1: 1.180041 0.960738 1.44940
    ##   2: 1.001086 0.892764 1.12255
    ##   3: 1.078753 0.931002 1.24995
    ##   4: 1.063509 0.939229 1.20423
    ##   5: 1.000378 0.915322 1.09334
    ##  ---                          
    ## 993: 0.974015 0.903328 1.05023
    ## 994: 1.025310 0.934141 1.12538
    ## 995: 0.948245 0.811418 1.10815
    ## 996: 1.001026 0.908670 1.10277
    ## 997: 0.966104 0.770456 1.21143

Acknowledgments
===============

*RegParallel* would not exist were it not for initial contributions from:

[Jessica Lasky-Su](https://connects.catalyst.harvard.edu/Profiles/display/Person/79780), [Myles Lewis](https://www.qmul.ac.uk/whri/people/academic-staff/items/lewismyles.html), [Michael Barnes](https://www.qmul.ac.uk/whri/people/academic-staff/items/barnesmichael.html)

Session info
============

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/atlas-base/atlas/libblas.so.3.0
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
    ##  [1] survminer_0.4.3             ggpubr_0.1.7               
    ##  [3] ggplot2_3.0.0               biomaRt_2.37.5             
    ##  [5] bindrcpp_0.2.2              GEOquery_2.49.0            
    ##  [7] stringr_1.3.1               DESeq2_1.21.10             
    ##  [9] magrittr_1.5                airway_0.115.0             
    ## [11] SummarizedExperiment_1.11.5 DelayedArray_0.7.35        
    ## [13] matrixStats_0.54.0          Biobase_2.41.2             
    ## [15] GenomicRanges_1.33.13       GenomeInfoDb_1.17.1        
    ## [17] IRanges_2.15.17             S4Vectors_0.19.19          
    ## [19] BiocGenerics_0.27.1         RegParallel_0.99.0         
    ## [21] survival_2.42-6             MASS_7.3-50                
    ## [23] data.table_1.11.6           doParallel_1.0.14          
    ## [25] doMC_1.3.5                  iterators_1.0.10           
    ## [27] foreach_1.4.4               BiocParallel_1.15.12       
    ## [29] knitr_1.20                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-137           cmprsk_2.2-7           bitops_1.0-6          
    ##  [4] bit64_0.9-7            httr_1.3.1             RColorBrewer_1.1-2    
    ##  [7] progress_1.2.0         rprojroot_1.3-2        tools_3.5.1           
    ## [10] backports_1.1.2        R6_2.2.2               rpart_4.1-13          
    ## [13] Hmisc_4.1-1            DBI_1.0.0              lazyeval_0.2.1        
    ## [16] colorspace_1.3-2       nnet_7.3-12            withr_2.1.2           
    ## [19] tidyselect_0.2.4       gridExtra_2.3          prettyunits_1.0.2     
    ## [22] bit_1.1-14             curl_3.2               compiler_3.5.1        
    ## [25] formatR_1.5            htmlTable_1.12         xml2_1.2.0            
    ## [28] labeling_0.3           scales_1.0.0           checkmate_1.8.5       
    ## [31] survMisc_0.5.5         readr_1.1.1            genefilter_1.63.0     
    ## [34] digest_0.6.16          foreign_0.8-70         rmarkdown_1.10        
    ## [37] XVector_0.21.3         base64enc_0.1-3        pkgconfig_2.0.1       
    ## [40] htmltools_0.3.6        highr_0.7              limma_3.37.3          
    ## [43] htmlwidgets_1.2        rlang_0.2.2            rstudioapi_0.7        
    ## [46] RSQLite_2.1.1          bindr_0.1.1            zoo_1.8-3             
    ## [49] acepack_1.4.1          dplyr_0.7.6            RCurl_1.95-4.11       
    ## [52] GenomeInfoDbData_1.1.0 Formula_1.2-3          Matrix_1.2-14         
    ## [55] Rcpp_0.12.18           munsell_0.5.0          stringi_1.2.4         
    ## [58] yaml_2.2.0             zlibbioc_1.27.0        plyr_1.8.4            
    ## [61] grid_3.5.1             blob_1.1.1             crayon_1.3.4          
    ## [64] lattice_0.20-35        splines_3.5.1          annotate_1.59.0       
    ## [67] hms_0.4.2              locfit_1.5-9.1         pillar_1.3.0          
    ## [70] geneplotter_1.59.0     codetools_0.2-15       XML_3.98-1.16         
    ## [73] glue_1.3.0             evaluate_0.11          latticeExtra_0.6-28   
    ## [76] gtable_0.2.0           purrr_0.2.5            tidyr_0.8.1           
    ## [79] km.ci_0.5-2            assertthat_0.2.0       broom_0.5.0           
    ## [82] xtable_1.8-2           tibble_1.4.2           KMsurv_0.1-5          
    ## [85] AnnotationDbi_1.43.1   memoise_1.1.0          cluster_2.0.7-1

References
----------

Blighe (2018)

Blighe, Kevin. 2018. “RegParallel: Standard regression functions in R enabled for parallel processing over large data-frames.” <https://github.com/kevinblighe>.
