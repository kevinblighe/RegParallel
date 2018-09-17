RegParallel <- function(
  data,
  formula,
  variables = NULL,
  blocksize = 1000,
  cores = 4,
  FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit'))
{
  blocksize <- round(blocksize, 0)
  message("Blocksize: ", blocksize)

  cores <- round(cores, 0)
  message("cores: ", cores * cores)
  require(BiocParallel)

  require(foreach)

  require(doMC)
  registerDoMC(cores)

  require(parallel)
  options("mc.cores"=cores)

  # UNIX
  require(doParallel)
  #registerDoParallel(makeCluster(detectCores() - 1))
  registerDoParallel(cores)

  require(data.table)
  
  # determine number of blocks
  blocks <- floor((length(variables)) / blocksize) + 1

  # create en empty list to hold the formulae
  formula.list <- list()

  # store each possible formula in the list
  for (i in 1:length(variables)) {
    formula.list[[i]] <- as.formula(paste(formula, " + ", variables[i], sep=""))
  }

  startIndex <- 9

  basemodel <- FUN(formula = as.formula(formula), data = data)

  foreach(l = 1:blocks, .combine = rbind, .multicombine = TRUE, .inorder = FALSE, .packages=c("data.table", "doParallel", "parallel", "doMC", "foreach", "BiocParallel")) %dopar% {
    # first block - will be executed just once
    if (l==1) {
      message(paste("Processing ", blocksize, " formulae, batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: 1; ", "Index2: ", (blocksize*l), sep=""))
      df <- data[,c(which(colnames(data) %in% c("dex", "cell")), startIndex + (1:(blocksize*l)))]
      models <- mclapply(formula.list[1:(blocksize*l)], function(f) summary(FUN(formula = f, data = df))$coefficients)
      names(models) <- variables[1:(blocksize*l)]

      # convert to data frames  
      wObjects <- mclapply(models, function(x) data.frame(rownames(x), x))

      # remove intercept and covariates from final output (OPTIONAL)
      wObjects <- mclapply(wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])

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

      return(wObjects)
    }

    # final block - will be executed just once
    if (l==blocks) {
      message(paste("Processing final batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: ", (1+(blocksize*(l-1))), "; ", "Index2: ", length(formula.list), sep=""))
      df <- data[,c(which(colnames(data) %in% c("dex", "cell")), startIndex + (((1+(blocksize*(l-1)))):(length(formula.list))))]
      models <- mclapply(formula.list[(1+(blocksize*(l-1))):length(formula.list)], function(f) summary(FUN(formula = f, data = df))$coefficients)
      names(models) <- variables[(1+(blocksize*(l-1))):length(formula.list)]

      # convert to data frames  
      wObjects <- mclapply(models, function(x) data.frame(rownames(x), x))

      # remove intercept and covariates from final output (OPTIONAL)
      wObjects <- mclapply(wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])

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

      return(wObjects)
    }

    # any other blocks - executed any number of times between first and final block
    if (l>1 && l<blocks) {
      message(paste("Processing ", blocksize, " formulae, batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: ", (1+(blocksize*(l-1))), "; ", "Index2: ", (blocksize*l), sep=""))
      df <- data[,c(which(colnames(data) %in% c("dex", "cell")), startIndex + ((1+(blocksize*(l-1))):(blocksize*l)))]
      models <- mclapply(formula.list[(1+(blocksize*(l-1))):(blocksize*l)], function(f) summary(FUN(formula = f, data = df))$coefficients)
      names(models) <- variables[(1+(blocksize*(l-1))):(blocksize*l)]

      # convert to data frames  
      wObjects <- mclapply(models, function(x) data.frame(rownames(x), x))

      # remove intercept and covariates from final output (OPTIONAL)
      wObjects <- mclapply(wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])

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

      return(wObjects)
    }
  }
}
