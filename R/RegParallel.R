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
  conflevel = 95)
{
  system <- Sys.info()['sysname']
  message('System is: ', system)

  blocksize <- round(blocksize, 0)
  message('Blocksize: ', blocksize)

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
    options('mc.cores' = cores)

    library(doMC, quietly = TRUE)
    registerDoMC(cores)
  }

  library(doParallel, quietly = TRUE)
  registerDoParallel(cores)

  library(data.table, quietly = TRUE)
  
  # determine number of blocks
  blocks <- floor((length(variables)) / blocksize) + 1

  # create en empty list to hold the formulae
  formula.list <- list()

  covs <- gsub(' ', '', unlist(strsplit(formula, '\\+|~')))
  covs <- covs[-grep('\\[x\\]', covs)]
  message('Variables included in model:')
  for (i in 1:length(covs)) {
    message('-- ', covs[i])
  }

  # store each possible formula in the list
  for (i in 1:length(variables)) {
    formula.list[[i]] <- as.formula(gsub('\\[x\\]', variables[i], formula))
  }

  five <- unlist(head(formula.list), 5)
  message("First 5 formulae:")
  for (i in 1:5) {
    message('-- ', five[i])
  }

  startIndex <- length(covs)

  data <- data[,which(colnames(data) %in% c(covs, variables))]

  foreach(l = 1:blocks, .combine = rbind, .multicombine = TRUE, .inorder = FALSE, .packages=c("data.table", "doParallel", "parallel", "doMC", "foreach", "BiocParallel")) %dopar% {
    # first block - will be executed just once
    if (l==1) {
      message(paste("Processing ", blocksize, " formulae, batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: 1; ", "Index2: ", (blocksize*l), sep=""))
      df <- data[,c(which(colnames(data) %in% c("dex", "cell")), startIndex + (1:(blocksize*l)))]

      if (nestedParallel == TRUE) {
        if (system == "Windows") {
            models <- parLapply(cl, formula.list[1:(blocksize*l)], function(f) summary(FUN(formula = f, data = df))$coefficients)
        } else {
            models <- mclapply(formula.list[1:(blocksize*l)], function(f) summary(FUN(formula = f, data = df))$coefficients)
        }
      } else if (nestedParallel == FALSE) {
        models <- lapply(formula.list[1:(blocksize*l)], function(f) summary(FUN(formula = f, data = df))$coefficients)
      } else {
        stop("Invalid value for argument nestedParallel. Must be TRUE/FALSE")
      }

      names(models) <- variables[1:(blocksize*l)]

      # convert to data frames
      if (system == "Windows") {
        wObjects <- parLapply(cl, models, function(x) data.frame(rownames(x), x))
      } else {
        wObjects <- mclapply(models, function(x) data.frame(rownames(x), x))
      }

      # remove intercept and covariates from final output
      if (system == "Windows") {
        wObjects <- parLapply(cl, wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])
      } else {
        wObjects <- mclapply(wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])
      }

      # detect failed models (will have 0 dimensions or return 'try' error)
      nullindices <- c(
        which(mapply(function(x) nrow(x)==0, wObjects)==TRUE),
        which(mapply(inherits, wObjects, 'try-error'))
      )
      for (i in nullindices) {
        wObjects[[i]] <- data.frame(t(c(names(wObjects)[i], NA, NA, NA, NA)))
      }

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), stringsAsFactors=FALSE)

      # set colnames
      colnames(wObjects) <- c("Variable", "Beta", "StandardError", "Z", "P")

      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(as.character(wObjects$StandardError))
      wObjects$Z <- as.numeric(as.character(wObjects$Z))
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate OR and confidence intervals
      wObjects$OR <- exp(wObjects[,'Beta'])
      wObjects$ORlower <- wObjects$OR * exp(qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$ORupperOR <- wObjects$OR * exp(qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      return(wObjects)
    }

    # final block - will be executed just once
    if (l==blocks) {
      message(paste("Processing final batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: ", (1+(blocksize*(l-1))), "; ", "Index2: ", length(formula.list), sep=""))
      df <- data[,c(which(colnames(data) %in% c("dex", "cell")), startIndex + (((1+(blocksize*(l-1)))):(length(formula.list))))]

      if (nestedParallel == TRUE) {
        if (system == "Windows") {
          models <- parLapply(cl, formula.list[(1+(blocksize*(l-1))):length(formula.list)], function(f) summary(FUN(formula = f, data = df))$coefficients)
        } else {
          models <- mclapply(formula.list[(1+(blocksize*(l-1))):length(formula.list)], function(f) summary(FUN(formula = f, data = df))$coefficients)
        }
      } else if (nestedParallel == FALSE) {
        models <- lapply(formula.list[(1+(blocksize*(l-1))):length(formula.list)], function(f) summary(FUN(formula = f, data = df))$coefficients)
      } else {
        stop("Invalid value for argument nestedParallel. Must be TRUE/FALSE")
      }

      names(models) <- variables[(1+(blocksize*(l-1))):length(formula.list)]

      # convert to data frames
      if (system == "Windows") {
        wObjects <- parLapply(cl, models, function(x) data.frame(rownames(x), x))
      } else {
        wObjects <- mclapply(models, function(x) data.frame(rownames(x), x))
      }

      # remove intercept and covariates from final output
      if (system == "Windows") {
        wObjects <- parLapply(cl, wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])
      } else {
        wObjects <- mclapply(wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])
      }

      # detect failed models (will have 0 dimensions or return 'try' error)
      nullindices <- c(
        which(mapply(function(x) nrow(x)==0, wObjects)==TRUE),
        which(mapply(inherits, wObjects, 'try-error'))
      )
      for (i in nullindices) {
        wObjects[[i]] <- data.frame(t(c(names(wObjects)[i], NA, NA, NA, NA)))
      }

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), stringsAsFactors=FALSE)

      # set colnames
      colnames(wObjects) <- c("Variable", "Beta", "StandardError", "Z", "P")

      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(as.character(wObjects$StandardError))
      wObjects$Z <- as.numeric(as.character(wObjects$Z))
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate OR and confidence intervals
      wObjects$OR <- exp(wObjects[,'Beta'])
      wObjects$ORlower <- wObjects$OR * exp(qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$ORupperOR <- wObjects$OR * exp(qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      return(wObjects)
    }

    # any other blocks - executed any number of times between first and final block
    if (l>1 && l<blocks) {
      message(paste("Processing ", blocksize, " formulae, batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: ", (1+(blocksize*(l-1))), "; ", "Index2: ", (blocksize*l), sep=""))
      df <- data[,c(which(colnames(data) %in% c("dex", "cell")), startIndex + ((1+(blocksize*(l-1))):(blocksize*l)))]

      if (nestedParallel == TRUE) {
        if (system == "Windows") {
          models <- parLapply(cl, formula.list[(1+(blocksize*(l-1))):(blocksize*l)], function(f) summary(FUN(formula = f, data = df))$coefficients)
        } else {
          models <- mclapply(formula.list[(1+(blocksize*(l-1))):(blocksize*l)], function(f) summary(FUN(formula = f, data = df))$coefficients)
        } 
      } else if (nestedParallel == FALSE) {
        models <- lapply(formula.list[(1+(blocksize*(l-1))):(blocksize*l)], function(f) summary(FUN(formula = f, data = df))$coefficients)
      } else {
        stop("Invalid value for argument nestedParallel. Must be TRUE/FALSE")
      }

      names(models) <- variables[(1+(blocksize*(l-1))):(blocksize*l)]

      # convert to data frames
      if (system == "Windows") {
        wObjects <- parLapply(cl, models, function(x) data.frame(rownames(x), x))
      } else {
        wObjects <- mclapply(models, function(x) data.frame(rownames(x), x))
      }

      # remove intercept and covariates from final output
      if (system == "Windows") {
        wObjects <- parLapply(cl, wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])
      } else {
        wObjects <- mclapply(wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])
      }

      # detect failed models (will have 0 dimensions or return 'try' error)
      nullindices <- c(
        which(mapply(function(x) nrow(x)==0, wObjects)==TRUE),
        which(mapply(inherits, wObjects, 'try-error'))
      )
      for (i in nullindices) {
        wObjects[[i]] <- data.frame(t(c(names(wObjects)[i], NA, NA, NA, NA)))
      }

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), stringsAsFactors=FALSE)

      # set colnames
      colnames(wObjects) <- c("Variable", "Beta", "StandardError", "Z", "P")

      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(as.character(wObjects$StandardError))
      wObjects$Z <- as.numeric(as.character(wObjects$Z))
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate OR and confidence intervals
      wObjects$OR <- exp(wObjects[,'Beta'])
      wObjects$ORlower <- wObjects$OR * exp(qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$ORupper <- wObjects$OR * exp(qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      return(wObjects)
    }
  }

  #if (system == 'Windows') {
  #  stopCluster(cl)
  #}
}
