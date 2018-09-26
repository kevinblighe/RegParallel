\name{coxphParallel}

\alias{coxphParallel}

\title{Standard regression functions in R enabled for parallel processing over large data-frames - Cox proportional hazards regression.}

\description{}

\usage{
coxphParallel(
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

  require(survival)
  data <- modelling[,1:800]
  variables <- colnames(data)[4:ncol(data)]
  data$time <- c(100,200,400,300,200,250,600,1000,886,450,
    c(100,200,400,300,200,250,600,1000,886,450)*1.5)
  data$alive <- c(0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,1,1,1,1)
  res4 <- RegParallel(
    data = data,
    formula = 'Surv(time, as.integer(alive)) ~ group * [*] + cell',
    FUN = function(formula, data)
      coxph(formula = formula,
        data = data,
        ties = 'breslow',
        singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = variables,
    blocksize = 399,
    cores = 2,
    nestedParallel = FALSE,
    conflevel = 97.5,
    excludeTerms = c('group', 'cell'),
    excludeIntercept = FALSE
  )

  # spot checks
  m <- coxph(formula = Surv(time, as.integer(alive)) ~ group * gene12 + cell, data = data, ties = 'breslow', singular.ok = TRUE)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.975)))
  res4[which(res4$Variable == 'gene12'),]

  m <- coxph(formula = Surv(time, as.integer(alive)) ~ group * gene267 + cell, data = data, ties = 'breslow', singular.ok = TRUE)
  summary(m)
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.975)))
  res4[which(res4$Variable == 'gene267'),]
}
