\name{clogitParallel}

\alias{clogitParallel}

\title{Standard regression functions in R enabled for parallel processing over large data-frames - conditional logistic regression.}

\description{}

\usage{
clogitParallel(
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
  res5 <- RegParallel(
    data = data,
    formula = 'as.integer(group) ~ [*] * strata(cell) + dosage',
    FUN = function(formula, data)
      clogit(formula = formula,
        data = data,
        ties = 'breslow',
        singular.ok = TRUE),
    FUNtype = 'clogit',
    variables = variables,
    blocksize = 200,
    cores = 2,
    nestedParallel = FALSE,
    conflevel = 50,
    excludeTerms = 'non-existent term',
    excludeIntercept = FALSE
  )

  # spot checks
  m <- clogit(formula = as.integer(group) ~ gene145 * strata(cell) + dosage, data = data, ties = 'breslow', singular.ok = TRUE)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.5)))
  res5[which(res5$Variable == 'gene145'),]

  m <- clogit(formula = as.integer(group) ~ gene34 * strata(cell) + dosage, data = data, ties = 'breslow', singular.ok = TRUE)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.5)))
  res5[which(res5$Variable == 'gene34'),]
}
