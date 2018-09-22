RegParallel <- function(
  data,
  formula,
  FUN = function(formula, data)
    glm(formula = formula,
      data = data,
      family = binomial(link = 'logit'),
      method = 'glm.fit'),
  FUNtype = 'glm',
  variables = NULL,
  blocksize = 1000,
  cores = 2,
  nestedParallel = FALSE,
  conflevel = 95,
  removeNULL = TRUE)
{
  system <- Sys.info()['sysname']
  message('System is: ', system)

  blocksize <- round(blocksize, 0)
  message('Blocksize: ', blocksize)

  if (blocksize > length(variables)) {
    stop(paste("blocksize is greater than number of variables to test.",
      "Choose a smaller blocksize"),
      call. = FALSE)
  }
  message('Cores / Threads: ', cores)

  if (nestedParallel == TRUE) {
    message('Nesting enabled. Total potential cores + forked threads is ', cores * cores)
  }

  library(BiocParallel, quietly = TRUE)

  library(foreach, quietly = TRUE)

  library(parallel, quietly = TRUE)

  if (system == 'Windows') {
    cl <- makeCluster(getOption('cl.cores', cores))
  } else {
    cl <- NULL
    options('mc.cores' = cores)

    library(doMC, quietly = TRUE)
    registerDoMC(cores)
  }

  library(doParallel, quietly = TRUE)
  registerDoParallel(cores)

  library(data.table, quietly = TRUE)

  # determine number of blocks
  blocks <- floor((length(variables)) / blocksize) + 1
  if (blocksize == length(variables)) {
    blocks <- 1
  }

  if (FUNtype == 'glm') {
    res <- glmParallel(
      data = data,
      formula = formula,
      FUN = FUN,
      variables = variables,
      blocksize = blocksize,
      blocks = blocks,
      system = system,
      nestedParallel = nestedParallel,
      conflevel = conflevel,
      removeNULL = removeNULL,
      cluster = cl)
  } else if (FUNtype == 'coxph') {
    res <- coxphParallel(
      data = data,
      formula = formula,
      FUN = FUN,
      variables = variables,
      blocksize = blocksize,
      blocks = blocks,
      system = system,
      nestedParallel = nestedParallel,
      conflevel = conflevel,
      removeNULL = removeNULL,
      cluster = cl)
  }

  return(res)
}
