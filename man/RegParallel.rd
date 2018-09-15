\name{RegParallel}

\alias{RegParallel}

\title{Parallelised regression functions}

\description{Parallelised regression functions}

\usage{
RegParallel(
  data,
  formula,
  variables,
  blocksize = 500,
  cores = 12,
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
  dds <- DESeqDataSet(airway, design = ~ cell + dex)
  dds <- DESeq(dds, betaPrior = FALSE)
  rlogcounts <- rlog(dds, blind=FALSE) # or vstcounts <- vst(dds, blind=TRUE)
  rlogcounts <- assay(rlogcounts)
  modelling <- data.frame(colData(airway), t(rlogcounts))

  rm(dds); rm(airway); rm(rlogcounts)

  data <- modelling[,1:500]

  res <- RegParallel(
    data = data,
    formula = "dex ~ cell",
    variables = colnames(data)[10:ncol(data)],
    blocksize = 200,
    cores = 3,
    FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit')
  )

