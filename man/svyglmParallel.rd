\name{svyglmParallel}

\alias{svyglmParallel}

\title{Standard regression functions in R enabled for parallel processing over large data-frames - generalised linear model, with survey weights}

\description{This is a non-user function that is managed by RegParallel, the primary function.}

\usage{
svyglmParallel(
  data,
  design,
  formula.list,
  FUN,
  variables,
  terms,
  startIndex,
  blocksize,
  blocks,
  APPLYFUN,
  conflevel,
  excludeTerms,
  excludeIntercept)
}

\arguments{
  \item{data}{A data-frame that contains all model terms to be tested.
  Variables that have all zeros will, automatically, be removed. REQUIRED.}
  \item{design}{A survey design, created by survey::svydesign. REQUIRED.}
  \item{formula.list}{A list containing formulae that can be coerced to
  formula class via as.formula(). REQUIRED.}
  \item{FUN}{Regression function. Must be of form, for example:
  function(formula, data) glm(formula = formula, family = binomial, data = data).
  REQUIRED.}
  \item{variables}{Vector of variable names in data to be tested
  independently. Each variable will have its own formula in formula.list.
  REQUIRED.}
  \item{terms}{Vector of terms used in the formulae in formula.list, excluding
  the primary variable of interest. REQUIRED.}
  \item{startIndex}{Starting column index in data object from which
  processing can commence. REQUIRED.}
  \item{blocksize}{Number of variables to test in each foreach loop.
  REQUIRED.}
  \item{blocks}{Total number of blocks required to complete analysis.
  REQUIRED.}
  \item{APPLYFUN}{The apply function to be used within each block during
  processing. Will be one of: 'mclapply(...)', system=linux/mac and
  nestedParallel=TRUE; 'parLapply(cl, ...)', system=windows and
  nestedParallel=TRUE; 'lapply(...)', nestedParallel=FALSE. REQUIRED.}
  \item{conflevel}{Confidence level for calculating odds or hazard ratios.
  REQUIRED.}
  \item{excludeTerms}{Remove these terms from the final output. These will
  simply be grepped out. REQUIRED.}
  \item{excludeIntercept}{Remove intercept terms from the final output.
  REQUIRED.}
}

\details{
This is a non-user function that is managed by RegParallel, the
primary function.
}

\value{
A \code{\link{data.table}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
}

\examples{
  require(survey)
  data(nhanes)
  design <- svydesign(id = ~ SDMVPSU,
    strata = ~ SDMVSTRA,
    weights = ~ WTMEC2YR,
    nest = TRUE,
    data = nhanes)
}
