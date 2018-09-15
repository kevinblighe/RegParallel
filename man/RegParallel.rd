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

data <- modelling

res <- RegParallel(
  data = data,
  formula = "dex ~ cell",
  variables = colnames(data)[10:ncol(data)],
  blocksize = 1000,
  cores = 26,
  FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit')
)

