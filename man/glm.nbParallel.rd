\name{glm.nbParallel}

\alias{glm.nbParallel}

\title{Standard regression functions in R enabled for parallel processing over large data-frames - negative binomial generalised linear model / negative binomial logistic regression.}

\description{}

\usage{
glm.nbParallel(
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

  data <- modelling[,1:5000]
  data[,4:ncol(data)] <- asinh(scale(data[,4:ncol(data)]))
  variables <- colnames(data)[4:ncol(data)]
  res7 <- RegParallel(
    data = data,
    formula = 'as.integer(cell) ~ [*] + group * dosage',
    FUN = function(formula, data)
      glm.nb(formula = formula,
        data = data),
    FUNtype = 'glm.nb',
    variables = variables,
    blocksize = 500,
    cores = 2,
    nestedParallel = FALSE,
    conflevel = 95,
    excludeTerms = NULL,
    excludeIntercept = FALSE
  )

  # spot checks
  m <- glm.nb(formula = as.integer(cell) ~ gene99 + group * dosage, data = data)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.95)))
  res7[which(res7$Variable == 'gene99'),]

  m <- glm.nb(formula = as.integer(cell) ~ gene2000 + group * dosage, data = data)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.95)))
  res7[which(res7$Variable == 'gene2000'),]
}
