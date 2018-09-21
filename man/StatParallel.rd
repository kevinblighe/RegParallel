\name{RegParallel}

\alias{RegParallel}

\title{Standard regression, correlation, and other functions in R enabled for parallel processing over large data-frames}

\description{In many analyses, a large amount of variables have to be tested independently against the trait/endpoint of interest, and also adjusted for covariates and confounding factors at the same time. The major botteleneck in these is the amount of time that it takes to complete these analyses.

With <i>StatParallel</i>, any number of tests can be performed simultaneously.  On a 12-core system, 144 variables can be tested simultaneously, with 1000s of variables processed in a matter of seconds.

Works for logistic regression, linear regression, conditional logistic regression, Cox proportional hazards and survival models, Bayesian logistic regression, ANOVA, and correlation analysis. Also works for GWAS studies loaded into R as snpMatrix objects.}

\usage{
RegParallel(
  data,
  formula,
  variables,
  blocksize = 1000,
  cores = 4,
  FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit'))
}

\arguments{
  \item{data}{A data-frame. REQUIRED.}
  \item{formula}{A formula. REQUIRED.}
  \item{variables}{Variables to test independently. REQUIRED.}
  \item{blocksize}{Number of variables to test in each foreach loop. DEFAULT = 500. OPTIONAL.}
  \item{cores}{CPU cores / threads. DEFAULT = 12. OPTIONAL.}
  \item{FUN}{Regression function. DEFAULT = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit'). OPTIONAL.}
}

\details{
  In many analyses, a large amount of variables have to be tested independently against the trait/endpoint of interest, and also adjusted for covariates and confounding factors at the same time. The major botteleneck in these is the amount of time that it takes to complete these analyses.

  With <i>StatParallel</i>, any number of tests can be performed simultaneously.  On a 12-core system, 144 variables can be tested simultaneously, with 1000s of variables processed in a matter of seconds.

  Works for logistic regression, linear regression, conditional logistic regression, Cox proportional hazards models, ANOVA, and correlations. Also works for GWAS studies loaded into R as snpMatrix objects.
}

\value{
A \code{\link{data-frame}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
}

\examples{

  library(airway)
  data("airway")
  library("DESeq2")
  dds <- DESeqDataSet(airway, design = ~ 1)
  dds <- DESeq(dds, betaPrior = FALSE)
  rlogcounts <- rlog(dds, blind=TRUE) # or vstcounts <- vst(dds, blind=TRUE)
  rlogcounts <- assay(rlogcounts)
  modelling <- data.frame(colData(airway), t(rlogcounts))

  rm(dds); rm(airway); rm(rlogcounts)

  data <- modelling[,1:4000]

  res <- RegParallel(
    data = data,
    formula = 'dex ~ [x]+cell',
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit'),
    FUNtype = 'glm',
    variables = colnames(data)[10:ncol(data)],
    blocksize = 700,
    cores = 2,
    nestedParallel = TRUE,
    conflevel = 99,
    removeNULL = TRUE
  )

  require(survival)
  data <- modelling[,1:3000]
  data$time <- c(12,6,13,5,12,7,21,5)
  data$alive <- c(1,0,1,0,1,0,1,0)
  res1 <- RegParallel(
    data = data,
    formula = 'Surv(alive) ~ [x] + strata(cell)',
    FUN = function(formula, data) coxph(formula = formula, data = data, ties = 'breslow', singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = colnames(modelling[,10:3000]),
    blocksize = 1000,
    cores = 2,
    nestedParallel = TRUE,
    removeNULL = FALSE
  )
}

