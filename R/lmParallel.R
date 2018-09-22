lmParallel <- function(
  data,
  formula,
  FUN,
  variables,
  blocksize,
  blocks,
  system,
  nestedParallel,
  conflevel,
  removeNULL,
  cluster)
{
  # create en empty list to hold the formulae
  formula.list <- list()

  # strip out information from the model to parse the variable names
  covs <- gsub(' ', '', unlist(strsplit(formula, '\\+|~|,')))

  # determine position of variable
  if (grep('\\[\\*\\]', covs) == 1) {
    xy <- 'y'
    message('Variable being tested is a dependent / y variable')
  } else {
    xy <- 'x'
    message('Variable being tested is a predictor / x variable')
  }

  covs <- covs[-grep('\\[\\*\\]', covs)]
  message('Variables included in model:')
  for (i in 1:length(covs)) {
    message('-- ', covs[i])
  }

  # store each possible formula in the list
  for (i in 1:length(variables)) {
    formula.list[[i]] <- as.formula(gsub('\\[\\*\\]', variables[i], formula))
  }

  five <- unlist(head(formula.list), 5)
  message("First 5 formulae:")
  for (i in 1:5) {
    message('-- ', five[i])
  }

  startIndex <- length(covs)

  # 'left align' the covariates
  data <- cbind(data[,which(colnames(data) %in% c(covs))], data[,which(colnames(data) %in% c(variables))])

  foreach(l = 1:blocks, .combine = rbind, .multicombine = TRUE, .inorder = FALSE, .packages = c("data.table", "doParallel", "parallel", "doMC", "foreach", "BiocParallel")) %dopar% {
    # first block - will be executed just once
    if (l==1) {
      message(paste("Processing ", blocksize, " formulae, batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: 1; ", "Index2: ", (blocksize*l), sep=""))

      df <- data[,c(which(colnames(data) %in% covs), startIndex + (1:(blocksize*l)))]

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
        wObjects <- parLapply(cl, names(models), function(x) data.frame(rep(x, length(rownames(models[[x]]))), rownames(models[[x]]), models[[x]], row.names=rownames(models[[x]])))
      } else {
        wObjects <- mclapply(names(models), function(x) data.frame(rep(x, length(rownames(models[[x]]))), rownames(models[[x]]), models[[x]], row.names=rownames(models[[x]])))
      }

      # remove intercept and covariates from final output
      if (system == "Windows") {
        if (xy == 'x') {
          wObjects <- parLapply(cl, wObjects, function(x) x[grep(paste(c("Intercept", covs), collapse="|^"), rownames(x), invert=TRUE),])
        } else if (xy == 'y') {
          wObjects <- parLapply(cl, wObjects, function(x) x[grep("Intercept", rownames(x), invert=TRUE),])
        }
      } else {
        if (xy == 'x') {
          wObjects <- mclapply(wObjects, function(x) x[grep(paste(c("Intercept", covs), collapse="|^"), rownames(x), invert=TRUE),])
        } else if (xy == 'y') {
          wObjects <- mclapply(wObjects, function(x) x[grep("Intercept", rownames(x), invert=TRUE),])
        }
      }

      # detect failed models (return 'try' error)
      nullindices <- c(
        which(mapply(inherits, wObjects, 'try-error'))
      )
      if ((removeNULL == TRUE) && length(nullindices) > 0) {
        wObjects <- wObjects[-nullindices]
      } else {
         for (i in nullindices) {
            wObjects[[i]] <- data.frame(t(c(names(wObjects)[i], rep(NA, cols))))
         }
      }

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), stringsAsFactors=FALSE)

      # set colnames
      colnames(wObjects) <- c('Variable', 'Term', 'Beta', 'StandardError', 't', 'P')

      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(as.character(wObjects$StandardError))
      wObjects$t <- as.numeric(as.character(wObjects$t))
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate OR and confidence intervals
      wObjects$OR <- exp(wObjects[,'Beta'])
      wObjects$ORlower <- wObjects$OR * exp(qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$ORupper <- wObjects$OR * exp(qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      # change Inf and NaN to NA
      wObjects$t[is.infinite(wObjects$t)] <- NA
      wObjects$t[wObjects$t == "NaN"] <- NA
      wObjects$P[is.infinite(wObjects$P)] <- NA
      wObjects$P[wObjects$P == "NaN"] <- NA
      wObjects$OR[is.infinite(wObjects$OR)] <- NA
      wObjects$OR[wObjects$OR == "NaN"] <- NA
      wObjects$ORlower[is.infinite(wObjects$ORlower)] <- NA
      wObjects$ORlower[wObjects$ORlower == "NaN"] <- NA
      wObjects$ORupper[is.infinite(wObjects$ORupper)] <- NA
      wObjects$ORupper[wObjects$ORupper == "NaN"] <- NA

      # if on Windows and there's only 1 block, then this is the first and final block
      # thus, we must free the cluster
      if ((blocks == 1) && (system == 'Windows')) {
        stopCluster(cl)
      }

      return(wObjects)
    }

    # final block - will be executed just once
    if (l==blocks) {
      message(paste("Processing final batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: ", (1+(blocksize*(l-1))), "; ", "Index2: ", length(formula.list), sep=""))
      df <- data[,c(which(colnames(data) %in% covs), startIndex + (((1+(blocksize*(l-1)))):(length(formula.list))))]

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
        wObjects <- parLapply(cl, names(models), function(x) data.frame(rep(x, length(rownames(models[[x]]))), rownames(models[[x]]), models[[x]], row.names=rownames(models[[x]])))
      } else {
        wObjects <- mclapply(names(models), function(x) data.frame(rep(x, length(rownames(models[[x]]))), rownames(models[[x]]), models[[x]], row.names=rownames(models[[x]])))
      }

      # remove intercept and covariates from final output
      if (system == "Windows") {
        if (xy == 'x') {
          wObjects <- parLapply(cl, wObjects, function(x) x[grep(paste(c("Intercept", covs), collapse="|^"), rownames(x), invert=TRUE),])
        } else if (xy == 'y') {
          wObjects <- parLapply(cl, wObjects, function(x) x[grep("Intercept", rownames(x), invert=TRUE),])
        }
      } else {
        if (xy == 'x') {
          wObjects <- mclapply(wObjects, function(x) x[grep(paste(c("Intercept", covs), collapse="|^"), rownames(x), invert=TRUE),])
        } else if (xy == 'y') {
          wObjects <- mclapply(wObjects, function(x) x[grep("Intercept", rownames(x), invert=TRUE),])
        }
      }

      # detect failed models (return 'try' error)
      nullindices <- c(
        which(mapply(inherits, wObjects, 'try-error'))
      )
      if ((removeNULL == TRUE) && length(nullindices) > 0) {
        wObjects <- wObjects[-nullindices]
      } else {
        for (i in nullindices) {
          wObjects[[i]] <- data.frame(t(c(names(wObjects)[i], rep(NA, cols))))
        }
      }

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), stringsAsFactors=FALSE)

      # set colnames
      colnames(wObjects) <- c('Variable', 'Term', 'Beta', 'StandardError', 't', 'P')

      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(as.character(wObjects$StandardError))
      wObjects$t <- as.numeric(as.character(wObjects$t))
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate OR and confidence intervals
      wObjects$OR <- exp(wObjects[,'Beta'])
      wObjects$ORlower <- wObjects$OR * exp(qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$ORupper <- wObjects$OR * exp(qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      # change Inf and NaN to NA
      wObjects$t[is.infinite(wObjects$t)] <- NA
      wObjects$t[wObjects$t == "NaN"] <- NA
      wObjects$P[is.infinite(wObjects$P)] <- NA
      wObjects$P[wObjects$P == "NaN"] <- NA
      wObjects$OR[is.infinite(wObjects$OR)] <- NA
      wObjects$OR[wObjects$OR == "NaN"] <- NA
      wObjects$ORlower[is.infinite(wObjects$ORlower)] <- NA
      wObjects$ORlower[wObjects$ORlower == "NaN"] <- NA
      wObjects$ORupper[is.infinite(wObjects$ORupper)] <- NA
      wObjects$ORupper[wObjects$ORupper == "NaN"] <- NA

      # final block. If Windows system, disable access to grabbed cluster
      if (system == 'Windows') {
        stopCluster(cl)
      }

      return(wObjects)
    }

    # any other blocks - executed any number of times between first and final block
    if (l>1 && l<blocks) {
      message(paste("Processing ", blocksize, " formulae, batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: ", (1+(blocksize*(l-1))), "; ", "Index2: ", (blocksize*l), sep=""))
      df <- data[,c(which(colnames(data) %in% covs), startIndex + ((1+(blocksize*(l-1))):(blocksize*l)))]

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
        wObjects <- parLapply(cl, names(models), function(x) data.frame(rep(x, length(rownames(models[[x]]))), rownames(models[[x]]), models[[x]], row.names=rownames(models[[x]])))
      } else {
        wObjects <- mclapply(names(models), function(x) data.frame(rep(x, length(rownames(models[[x]]))), rownames(models[[x]]), models[[x]], row.names=rownames(models[[x]])))
      }

      # remove intercept and covariates from final output
      if (system == "Windows") {
        if (xy == 'x') {
          wObjects <- parLapply(cl, wObjects, function(x) x[grep(paste(c("Intercept", covs), collapse="|^"), rownames(x), invert=TRUE),])
        } else if (xy == 'y') {
          wObjects <- parLapply(cl, wObjects, function(x) x[grep("Intercept", rownames(x), invert=TRUE),])
        }
      } else {
        if (xy == 'x') {
          wObjects <- mclapply(wObjects, function(x) x[grep(paste(c("Intercept", covs), collapse="|^"), rownames(x), invert=TRUE),])
        } else if (xy == 'y') {
          wObjects <- mclapply(wObjects, function(x) x[grep("Intercept", rownames(x), invert=TRUE),])
        }
      }

      # detect failed models (return 'try' error)
      nullindices <- c(
        which(mapply(inherits, wObjects, 'try-error'))
      )
      if ((removeNULL == TRUE) && length(nullindices) > 0) {
        wObjects <- wObjects[-nullindices]
      } else {
        for (i in nullindices) {
          wObjects[[i]] <- data.frame(t(c(names(wObjects)[i], rep(NA, cols))))
        }
      }

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), stringsAsFactors=FALSE)

      # set colnames
      colnames(wObjects) <- c('Variable', 'Term', 'Beta', 'StandardError', 't', 'P')

      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(as.character(wObjects$StandardError))
      wObjects$t <- as.numeric(as.character(wObjects$t))
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate OR and confidence intervals
      wObjects$OR <- exp(wObjects[,'Beta'])
      wObjects$ORlower <- wObjects$OR * exp(qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$ORupper <- wObjects$OR * exp(qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      # change Inf and NaN to NA
      wObjects$t[is.infinite(wObjects$t)] <- NA
      wObjects$t[wObjects$t == "NaN"] <- NA
      wObjects$P[is.infinite(wObjects$P)] <- NA
      wObjects$P[wObjects$P == "NaN"] <- NA
      wObjects$OR[is.infinite(wObjects$OR)] <- NA
      wObjects$OR[wObjects$OR == "NaN"] <- NA
      wObjects$ORlower[is.infinite(wObjects$ORlower)] <- NA
      wObjects$ORlower[wObjects$ORlower == "NaN"] <- NA
      wObjects$ORupper[is.infinite(wObjects$ORupper)] <- NA
      wObjects$ORupper[wObjects$ORupper == "NaN"] <- NA

      return(wObjects)
    }
  }
}
