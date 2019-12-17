\name{RegParallel}

\alias{RegParallel}

\title{Standard regression functions in R enabled for parallel processing over large data-frames}

\description{In many analyses, a large amount of variables have to be tested independently against the trait/endpoint of interest, and also adjusted for covariates and confounding factors at the same time. The major botteleneck in these is the amount of time that it takes to complete these analyses.

With RegParallel, a large number of tests can be performed simultaneously.  On a 12-core system, 144 variables can be tested simultaneously, with 1000s of variables processed in a matter of seconds via 'nested' parallel processing.

Works for logistic regression, linear regression, conditional logistic regression, Cox proportional hazards and survival models, Bayesian logistic regression, and negative binomial regression.}

\usage{
RegParallel(
  data,
  formula,
  FUN,
  FUNtype,
  variables,
  blocksize = 500,
  cores = 3,
  nestedParallel = FALSE,
  p.adjust = 'none',
  conflevel = 95,
  excludeTerms = NULL,
  excludeIntercept = TRUE)
}

\arguments{
  \item{data}{A data-frame that contains all model terms to be tested.
  Variables that have all zeros will, automatically, be removed. REQUIRED.}
  \item{formula}{A valid formula. Excluding the '[*]' term, which is reserved
  for RegParallel and indicates the position in the formula for the variable
  of interest, must pass as.formula() check. REQUIRED.}
  \item{FUN}{Regression function. Must be of form, for example:
  function(formula, data) glm(formula = formula, family=binomial, data = data).
  REQUIRED.}
  \item{FUNtype}{Regression function type. Must be one of 'glm', 'lm',
  'coxph', 'clogit', 'bayesglm', or 'glm.nb'. REQUIRED.}
  \item{variables}{Vector of variable names in data to be tested
  independently. Each variable will take the place of '[*]' in the supplied
  formula. REQUIRED.}
  \item{blocksize}{Number of variables to test in each foreach loop.
  DEFAULT = 500. OPTIONAL.}
  \item{cores}{CPU cores / threads. DEFAULT = 3. OPTIONAL.}
  \item{nestedParallel}{In RegParallel, parallelisation initially occurs at
  the block level, ie., multiple blocks of models are processed in parallel.
  If nestedParallel is enabled, a second level of parallelisation occurs
  within each block in addition. Warning! - this doubles the usage of cores.
  DEFAULT = FALSE. OPTIONAL.}
  \item{p.adjust}{Method for adjusting p-values for false discovery rate. Must
  be one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr',
  'none'. See ?p.adjust for further details. DEFAULT = 'none'. OPTIONAL}
  \item{conflevel}{Confidence level for calculating odds or hazard ratios.
  DEFAULT = 95. OPTIONAL.}
  \item{excludeTerms}{Remove these terms from the final output. These will
  simply be grepped out. DEFAULT = NULL. OPTIONAL.}
  \item{excludeIntercept}{Remove intercept terms from the final output.
  DEFAULT = TRUE. OPTIONAL.}
}

\details{
In many analyses, a large amount of variables have to be tested independently against the trait/endpoint of interest, and also adjusted for covariates and confounding factors at the same time. The major bottleneck in these is the amount of time that it takes to complete these analyses. With RegParallel, a large number of tests can be performed simultaneously. On a 12-core system, 144 variables can be tested simultaneously, with 1000s of variables processed in a matter of seconds via 'nested' parallel processing. Works for logistic regression, linear regression, conditional logistic regression, Cox proportional hazards and survival models, and Bayesian logistic regression.
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
    #p.adjust = "bonferroni",
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
    p.adjust = "none",
    conflevel = 90,
    excludeTerms = NULL,
    excludeIntercept = FALSE
  )

  # spot checks
  m <- glm(gene29 ~ cell:dosage, data=data, family=gaussian)
  summary(m)$coefficients
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.90)))
  res2[which(res2$Variable == 'gene29'),]


  ###


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
    p.adjust = "holm",
    conflevel = 99.999,
    excludeTerms = NULL,
    excludeIntercept = FALSE
  )

  # spot checks
  m <- lm(as.numeric(gene454) ~ dosage ^ 3, data=data)
  summary(m)$coefficients
  exp(cbind("Odds ratio" = coef(m), confint.default(m, level = 0.99999)))
  res3[which(res3$Variable == 'gene454'),]


  ###


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
    p.adjust = "hommel",
    conflevel = 97.5,
    excludeTerms = c('group', 'cell'),
    excludeIntercept = FALSE
  )

  # spot checks
  m <- coxph(formula = Surv(time, as.integer(alive)) ~ group * gene12 + cell, data = data, ties = 'breslow', singular.ok = TRUE)
  summary(m)
  exp(cbind("Hazards ratio" = coef(m), confint.default(m, level = 0.975)))
  res4[which(res4$Variable == 'gene12'),]

  m <- coxph(formula = Surv(time, as.integer(alive)) ~ group * gene267 + cell, data = data, ties = 'breslow', singular.ok = TRUE)
  summary(m)
  exp(cbind("Hazards ratio" = coef(m), confint.default(m, level = 0.975)))
  res4[which(res4$Variable == 'gene267'),]


  ###


  data <- modelling[,1:500]
  variables <- colnames(data)[4:ncol(data)]
  res5 <- RegParallel(
    data = data,
    formula = 'as.integer(factor(group)) ~ [*] * strata(cell) + dosage',
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
    p.adjust = "fdr",
    conflevel = 50,
    excludeTerms = 'non-existent term',
    excludeIntercept = FALSE
  )

  # spot checks
  m <- clogit(formula = as.integer(factor(group)) ~ gene145 * strata(cell) + dosage, data = data, ties = 'breslow', singular.ok = TRUE)
  summary(m)
  exp(cbind("Hazards ratio" = coef(m), confint.default(m, level = 0.5)))
  res5[which(res5$Variable == 'gene145'),]

  m <- clogit(formula = as.integer(factor(group)) ~ gene34 * strata(cell) + dosage, data = data, ties = 'breslow', singular.ok = TRUE)
  summary(m)
  exp(cbind("Hazards ratio" = coef(m), confint.default(m, level = 0.5)))
  res5[which(res5$Variable == 'gene34'),]


  ###


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
    p.adjust = "fdr",
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
