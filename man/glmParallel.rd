\name{glmParallel}

\alias{glmParallel}

\title{Standard regression functions in R enabled for parallel processing over large data-frames - generalised linear model / logistic regression.}

\description{}

\usage{
glmParallel(
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

  data <- modelling[,1:2000]
  variables <- colnames(data)[4:ncol(data)]
  res1 <- RegParallel(
    data = data,
    formula = 'factor(group) ~ [*] + (cell:dosage) ^ 2',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = binomial(link = 'logit'),
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 700,
    cores = 2,
    nestedParallel = TRUE,
    conflevel = 99,
    excludeTerms = NULL,
    excludeIntercept = TRUE
  )

  # spot checks
  m <- glm(factor(group) ~ gene265 + (cell:dosage) ^ 2, data=data, family=binomial)
  summary(m)$coefficients
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.99)))
  res1[which(res1$Variable == 'gene265'),]

  m <- glm(factor(group) ~ gene1688 + (cell:dosage) ^ 2, data=data, family=binomial)
  summary(m)$coefficients
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.99)))
  res1[which(res1$Variable == 'gene1688'),]


  ###


  data <- modelling[,1:500]
  variables <- colnames(data)[4:ncol(data)]
  res2 <- RegParallel(
    data = data,
    formula = '[*] ~ cell:dosage',
    FUN = function(formula, data)
      glm(formula = formula,
        data = data,
        family = gaussian,
        method = 'glm.fit'),
    FUNtype = 'glm',
    variables = variables,
    blocksize = 496,
    cores = 2,
    nestedParallel = TRUE,
    conflevel = 90,
    excludeTerms = NULL,
    excludeIntercept = FALSE
  )

  # spot checks
  m <- glm(gene29 ~ cell:dosage, data=data, family=gaussian)
  summary(m)$coefficients
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.90)))
  res2[which(res2$Variable == 'gene29'),]
}
