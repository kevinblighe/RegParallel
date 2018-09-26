\name{bayesglmParallel}

\alias{bayesglmParallel}

\title{Standard regression functions in R enabled for parallel processing over large data-frames - Bayesian logistic regression}

\description{}

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
  system,
  cluster,
  nestedParallel,
  conflevel,
  excludeTerms,
  excludeIntercept)
}

\arguments{}

\details{}

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
    formula = 'as.numeric(cell) ~ [*]:dosage',
    FUN = function(formula, data)
      bayesglm(formula = formula,
        data = data,
        prior.mean = 2),
    FUNtype = 'bayesglm',
    variables = variables,
    blocksize = 500,
    cores = 2,
    nestedParallel = FALSE,
    conflevel = 99,
    excludeTerms = NULL,
    excludeIntercept = FALSE
  )

  # spot checks
  m <- bayesglm(formula = as.numeric(cell) ~ gene1645:dosage, data = data, prior.mean = 2)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.99)))
  res6[which(res6$Variable == 'gene1645'),]

  m <- bayesglm(formula = as.numeric(cell) ~ gene3664:dosage, data = data, prior.mean = 2)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.99)))
  res6[which(res6$Variable == 'gene3664'),]
}
