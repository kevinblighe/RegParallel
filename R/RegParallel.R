RegParallel <- function(
  data,
  formula,
  variables = NULL,
  blocksize = 500,
  cores = 2,
  FUN = function(formula, data) glm(formula = formula, data = data, family = binomial(link = 'logit'), method = 'glm.fit'))
{

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
      models <- mclapply(formula.list[1:(blocksize*l)], function(f) FUN(formula = f, data = df))
      names(models) <- variables[1:(blocksize*l)]

      #If any models failed, detect them and replace with the base model
      wObjects <- mclapply(models, function(x) if (is.null(names(x))) { x <- basemodel } else {x} )

      #Extract information from the model and store in a new list
      wObjects <- mclapply(wObjects, function(x) data.frame(rownames(summary(x)$coefficients), summary(x)$coefficients))

      # remove intercept and covariates from final output (OPTIONAL)
      wObjects <- mclapply(wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])

      #Convert the list into a data frame for writing
      wObject <- do.call(rbind, mclapply(wObjects, data.frame, stringsAsFactors=FALSE))
      wObject <- data.frame(gsub("\\.[a-zA-Z0-9_():]*", "", rownames(wObject)), wObject)
      wObject[,2] <- gsub("factor", "", wObject[,2])

      return(data.table(wObject[,-c(1)]))
    }

    # final block - will be executed just once
    if (l==blocks) {
      message(paste("Processing final batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: ", (1+(blocksize*(l-1))), "; ", "Index2: ", length(formula.list), sep=""))
      df <- data[,c(which(colnames(data) %in% c("dex", "cell")), startIndex + (((1+(blocksize*(l-1)))):(length(formula.list))))]
      models <- mclapply(formula.list[(1+(blocksize*(l-1))):length(formula.list)], function(x) glm(x, df, family=binomial(link='logit')))
      names(models) <- variables[(1+(blocksize*(l-1))):length(formula.list)]

      #If any models failed, detect them and replace with the base model
      wObjects <- mclapply(models, function(x) if (is.null(names(x))) { x <- basemodel } else {x} )

      #Extract information from the model and store in a new list
      wObjects <- mclapply(wObjects, function(x) data.frame(rownames(summary(x)$coefficients), summary(x)$coefficients))

      # remove intercept and covariates from final output (OPTIONAL)
      wObjects <- mclapply(wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])

      #Convert the list into a data frame for writing
      wObject <- do.call(rbind, mclapply(wObjects, data.frame, stringsAsFactors=FALSE))
      wObject <- data.frame(gsub("\\.[a-zA-Z0-9_():]*", "", rownames(wObject)), wObject)
      wObject[,2] <- gsub("factor", "", wObject[,2])

      return(data.table(wObject[,-c(1)]))
    }

    # any other blocks - executed any number of times between first and final block
    if (l>1 && l<blocks) {
      message(paste("Processing ", blocksize, " formulae, batch ", l, " of ", blocks, sep=""))
      message(paste("Index1: ", (1+(blocksize*(l-1))), "; ", "Index2: ", (blocksize*l), sep=""))
      df <- data[,c(which(colnames(data) %in% c("dex", "cell")), startIndex + ((1+(blocksize*(l-1))):(blocksize*l)))]
      models <- mclapply(formula.list[(1+(blocksize*(l-1))):(blocksize*l)], function(f) FUN(formula = f, data = df))
      names(models) <- variables[(1+(blocksize*(l-1))):(blocksize*l)]

      #If any models failed, detect them and replace with the base model
      wObjects <- mclapply(models, function(x) if (is.null(names(x))) { x <- basemodel } else {x} )

      #Extract information from the model and store in a new list
      wObjects <- mclapply(wObjects, function(x) data.frame(rownames(summary(x)$coefficients), summary(x)$coefficients))

      # remove intercept and covariates from final output (OPTIONAL)
      wObjects <- mclapply(wObjects, function(x) x[grep("Intercept|^cell", rownames(x), invert=TRUE),])

      #Convert the list into a data frame for writing
      wObject <- do.call(rbind, mclapply(wObjects, data.frame, stringsAsFactors=FALSE))
      wObject <- data.frame(gsub("\\.[a-zA-Z0-9_():]*", "", rownames(wObject)), wObject)
      wObject[,2] <- gsub("factor", "", wObject[,2])

      return(data.table(wObject[,-c(1)]))
    }
  }
}
