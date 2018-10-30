\name{bayesglmParallel}

\alias{bayesglmParallel}

\title{Standard regression functions in R enabled for parallel processing over large data-frames - Bayesian logistic regression}

\description{This is a non-user function that is managed by RegParallel, the primary function.}

\usage{
bayesglmParallel(
  data,
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

  options(scipen=10)
  options(digits=6)

  col <- 20000
  row <- 20
  mat <- matrix(
    rexp(col*row, rate = .1),
    ncol = col)
  colnames(mat) <- paste0('gene', 1:ncol(mat))
  rownames(mat) <- paste0('sample', 1:nrow(mat))

  modelling <- data.frame(
    cell = rep(c('B', 'T'), nrow(mat) / 2),
    group = c(rep(c('treatment'), nrow(mat) / 2), rep(c('control'), nrow(mat) / 2)),
    dosage = t(data.frame(matrix(rexp(row, rate = 1), ncol = row))),
    mat,
    row.names = rownames(mat))

  data <- modelling[,1:5000]
  variables <- colnames(data)[4:ncol(data)]
  res6 <- RegParallel(
    data = data,
    formula = 'as.numeric(factor(cell)) ~ [*]:dosage',
    FUN = function(formula, data)
      bayesglm(formula = formula,
        data = data,
        prior.mean = 2),
    FUNtype = 'bayesglm',
    variables = variables,
    blocksize = 500,
    cores = 2,
    nestedParallel = FALSE,
    p.adjust = "none",
    conflevel = 99,
    excludeTerms = NULL,
    excludeIntercept = FALSE
  )

  # spot checks
  m <- bayesglm(formula = as.numeric(factor(cell)) ~ gene1645:dosage, data = data, prior.mean = 2)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.99)))
  res6[which(res6$Variable == 'gene1645'),]

  m <- bayesglm(formula = as.numeric(factor(cell)) ~ gene3664:dosage, data = data, prior.mean = 2)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.99)))
  res6[which(res6$Variable == 'gene3664'),]
}
