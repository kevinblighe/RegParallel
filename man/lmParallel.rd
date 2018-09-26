\name{lmParallel}

\alias{lmParallel}

\title{Standard regression functions in R enabled for parallel processing over large data-frames - linear model.}

\description{}

\usage{
lmParallel(
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
  excludeTerms)
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

  data <- modelling[,1:500]
  variables <- colnames(data)[4:ncol(data)]
  res3 <- RegParallel(
    data = data,
    formula = 'as.numeric([*]) ~ dosage ^ 3',
    FUN = function(formula, data)
      lm(formula = formula,
        data = data),
    FUNtype = 'lm',
    variables = variables,
    blocksize = 200,
    cores = 2,
    nestedParallel = FALSE,
    conflevel = 99.999,
    excludeTerms = NULL,
    excludeIntercept = FALSE
  )

  # spot checks
  m <- lm(as.numeric(gene454) ~ dosage ^ 3, data=data)
  summary(m)$coefficients
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.99999)))
  res3[which(res3$Variable == 'gene454'),]
}
