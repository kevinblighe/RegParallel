glmParallel <- function(
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
{

  ExpBeta <- l <- NULL

  foreach(l = 1:blocks,
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table',
      'doParallel', 'parallel',
      'foreach', 'BiocParallel')) %dopar% {

    # first block - will be executed just once
    if (l==1) {
      message(paste('Processing ', blocksize,
        ' formulae, batch ', l, ' of ', blocks, sep=''))
      message(paste('-- index1: 1; ', 'index2: ', (blocksize*l), sep=''))

      df <- data[,c(
        which(colnames(data) %in% terms),
        startIndex + (1:(blocksize*l)))]

      if (nestedParallel == TRUE) {
        if (system == 'Windows') {
            models <- parLapply(cluster, formula.list[1:(blocksize*l)],
              function(f) summary(FUN(formula = f, data = df))$coefficients)
        } else {
            models <- mclapply(formula.list[1:(blocksize*l)],
              function(f) summary(FUN(formula = f, data = df))$coefficients)
        }
      } else if (nestedParallel == FALSE) {
        models <- lapply(formula.list[1:(blocksize*l)],
          function(f) summary(FUN(formula = f, data = df))$coefficients)
      } else {
        stop('Invalid value for argument nestedParallel. Must be TRUE/FALSE')
      }

      names(models) <- variables[1:(blocksize*l)]

      # detect failed models (return 'try' error)
      nullindices <- c(
        which(mapply(inherits, models, 'try-error'))
      )
      if (length(nullindices) > 1) {
        print(models[[nullindices[1]]][1])
        stop(paste('\nOne or more models have failed. ',
          'Please recheck your data and model formula. ',
          'In particular, check that your variable names ',
          'are valid, i.e., no spaces, do not begin with ',
          'a number, no hyphens or other punctuation, et ',
          'cetera.', sep=''))
      }

      # convert to data frames
      if (system == 'Windows') {
        wObjects <- parLapply(cluster, names(models),
          function(x) data.frame(
            rep(x, length(rownames(models[[x]]))),
            rownames(models[[x]]),
            models[[x]],
            row.names=rownames(models[[x]])))
      } else {
        wObjects <- mclapply(names(models),
          function(x) data.frame(
            rep(x, length(rownames(models[[x]]))),
            rownames(models[[x]]),
            models[[x]],
            row.names=rownames(models[[x]])))
      }

      # remove intercept and / or specified terms from the output
      if (!is.null(excludeTerms) && excludeIntercept == TRUE) {
        if (system == 'Windows') {
          wObjects <- parLapply(cluster, wObjects,
            function(x) x[grep(
              paste(c('\\(Intercept\\)', excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        } else {
          wObjects <- mclapply(wObjects,
            function(x) x[grep(
              paste(c('\\(Intercept\\)', excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        }
      } else if (is.null(excludeTerms) && excludeIntercept == TRUE) {
        if (system == 'Windows') {
          wObjects <- parLapply(cluster, wObjects,
            function(x) x[grep(
              '\\(Intercept\\)',
              rownames(x),
              invert=TRUE),])
        } else {
          wObjects <- mclapply(wObjects,
            function(x) x[grep(
              '\\(Intercept\\)',
              rownames(x),
              invert=TRUE),])
        }
      } else if (!is.null(excludeTerms) && excludeIntercept == FALSE) {
        if (system == 'Windows') {
          wObjects <- parLapply(cluster, wObjects,
            function(x) x[grep(
              paste(c(excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        } else {
          wObjects <- mclapply(wObjects,
            function(x) x[grep(
              paste(c(excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        }
      }

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), stringsAsFactors=FALSE)

      # set colnames
      # detect if the y variable in a glm is a continous variable
      # (model will contain t stat; otherwise, Z)
      if (any(grepl('^t\\.value$', colnames(wObjects))) == TRUE) {
        colnames(wObjects) <- c('Variable', 'Term', 'Beta',
          'StandardError', 't', 'P')
      } else {
        colnames(wObjects) <- c('Variable', 'Term', 'Beta',
          'StandardError', 'Z', 'P')
      }

      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(
        as.character(wObjects$StandardError))
      if (any(grepl('^t$', colnames(wObjects))) == TRUE) {
        wObjects$t <- as.numeric(as.character(wObjects$t))
      } else {
        wObjects$Z <- as.numeric(as.character(wObjects$Z))
      }
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate OR and confidence intervals
      wObjects$OR <- exp(wObjects[,'Beta'])
      wObjects$ORlower <- wObjects$OR * exp(
        qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$ORupper <- wObjects$OR * exp(
        qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      # change Inf and NaN to NA
      if (any(grepl('^t$', colnames(wObjects)) == TRUE)) {
        wObjects$t[is.infinite(wObjects$t)] <- NA
        wObjects$t[wObjects$t == 'NaN'] <- NA
      } else {
        wObjects$Z[is.infinite(wObjects$Z)] <- NA
        wObjects$Z[wObjects$Z == 'NaN'] <- NA
      }
      wObjects$P[is.infinite(wObjects$P)] <- NA
      wObjects$P[wObjects$P == 'NaN'] <- NA
      wObjects$OR[is.infinite(wObjects$OR)] <- NA
      wObjects$OR[wObjects$OR == 'NaN'] <- NA
      wObjects$ORlower[is.infinite(wObjects$ORlower)] <- NA
      wObjects$ORlower[wObjects$ORlower == 'NaN'] <- NA
      wObjects$ORupper[is.infinite(wObjects$ORupper)] <- NA
      wObjects$ORupper[wObjects$ORupper == 'NaN'] <- NA

      # if on Windows and there's only 1 block,
      # then this is the first and final block;
      # thus, we must free the cluster
      if ((blocks == 1) && (system == 'Windows')) {
        stopCluster(cl)
      }

      return(wObjects)
    }

    # final block - will be executed just once
    if (l==blocks) {
      message(paste('Processing final batch ', l, ' of ',
        blocks, sep=''))
      message(paste('-- index1: ', (1+(blocksize*(l-1))), '; ',
        'index2: ', length(formula.list), sep=''))
      df <- data[,c(
        which(colnames(data) %in% terms),
        startIndex + (((1+(blocksize*(l-1)))):(length(formula.list))))]

      if (nestedParallel == TRUE) {
        if (system == 'Windows') {
          models <- parLapply(cl,
            formula.list[(1+(blocksize*(l-1))):length(formula.list)],
            function(f) summary(FUN(formula = f, data = df))$coefficients)
        } else {
          models <- mclapply(
            formula.list[(1+(blocksize*(l-1))):length(formula.list)],
            function(f) summary(FUN(formula = f, data = df))$coefficients)
        }
      } else if (nestedParallel == FALSE) {
        models <- lapply(
          formula.list[(1+(blocksize*(l-1))):length(formula.list)],
          function(f) summary(FUN(formula = f, data = df))$coefficients)
      } else {
        stop('Invalid value for argument nestedParallel. Must be TRUE/FALSE')
      }

      names(models) <- variables[(1+(blocksize*(l-1))):length(formula.list)]

      # detect failed models (return 'try' error)
      nullindices <- c(
        which(mapply(inherits, models, 'try-error'))
      )
      if (length(nullindices) > 1) {
        print(models[[nullindices[1]]][1])
        stop(paste('\nOne or more models have failed. ',
          'Please recheck your data and model formula. ',
          'In particular, check that your variable names ',
          'are valid, i.e., no spaces, do not begin with ',
          'a number, no hyphens or other punctuation, et ',
          'cetera.', sep=''))
      }

      # convert to data frames
      if (system == 'Windows') {
        wObjects <- parLapply(cluster, names(models),
          function(x) data.frame(
            rep(x, length(rownames(models[[x]]))),
            rownames(models[[x]]),
            models[[x]],
            row.names=rownames(models[[x]])))
      } else {
        wObjects <- mclapply(names(models),
          function(x) data.frame(
            rep(x, length(rownames(models[[x]]))),
            rownames(models[[x]]),
            models[[x]],
            row.names=rownames(models[[x]])))
      }

      # remove intercept and / or specified terms from the output
      if (!is.null(excludeTerms) && excludeIntercept == TRUE) {
        if (system == 'Windows') {
          wObjects <- parLapply(cluster, wObjects,
            function(x) x[grep(
              paste(c('\\(Intercept\\)', excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        } else {
          wObjects <- mclapply(wObjects,
            function(x) x[grep(
              paste(c('\\(Intercept\\)', excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        }
      } else if (is.null(excludeTerms) && excludeIntercept == TRUE) {
        if (system == 'Windows') {
          wObjects <- parLapply(cluster, wObjects,
            function(x) x[grep(
              '\\(Intercept\\)',
              rownames(x),
              invert=TRUE),])
        } else {
          wObjects <- mclapply(wObjects,
            function(x) x[grep(
              '\\(Intercept\\)',
              rownames(x),
              invert=TRUE),])
        }
      } else if (!is.null(excludeTerms) && excludeIntercept == FALSE) {
        if (system == 'Windows') {
          wObjects <- parLapply(cluster, wObjects,
            function(x) x[grep(
              paste(c(excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        } else {
          wObjects <- mclapply(wObjects,
            function(x) x[grep(
              paste(c(excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        }
      }

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), stringsAsFactors=FALSE)

      # set colnames
      # detect if the y variable in a glm is a continous variable
      # (model will contain t stat; otherwise, Z)
      if (any(grepl('^t\\.value$', colnames(wObjects))) == TRUE) {
        colnames(wObjects) <- c('Variable', 'Term', 'Beta',
          'StandardError', 't', 'P')
      } else {
        colnames(wObjects) <- c('Variable', 'Term', 'Beta',
          'StandardError', 'Z', 'P')
      }

      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(
        as.character(wObjects$StandardError))
      if (any(grepl('^t$', colnames(wObjects))) == TRUE) {
        wObjects$t <- as.numeric(as.character(wObjects$t))
      } else {
        wObjects$Z <- as.numeric(as.character(wObjects$Z))
      }
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate OR and confidence intervals
      wObjects$OR <- exp(wObjects[,'Beta'])
      wObjects$ORlower <- wObjects$OR * exp(
        qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$ORupper <- wObjects$OR * exp(
        qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      # change Inf and NaN to NA
      if (any(grepl('^t$', colnames(wObjects)) == TRUE)) {
        wObjects$t[is.infinite(wObjects$t)] <- NA
        wObjects$t[wObjects$t == 'NaN'] <- NA
      } else {
        wObjects$Z[is.infinite(wObjects$Z)] <- NA
        wObjects$Z[wObjects$Z == 'NaN'] <- NA
      }
      wObjects$P[is.infinite(wObjects$P)] <- NA
      wObjects$P[wObjects$P == 'NaN'] <- NA
      wObjects$OR[is.infinite(wObjects$OR)] <- NA
      wObjects$OR[wObjects$OR == 'NaN'] <- NA
      wObjects$ORlower[is.infinite(wObjects$ORlower)] <- NA
      wObjects$ORlower[wObjects$ORlower == 'NaN'] <- NA
      wObjects$ORupper[is.infinite(wObjects$ORupper)] <- NA
      wObjects$ORupper[wObjects$ORupper == 'NaN'] <- NA

      # final block. If Windows system, disable access to grabbed cluster
      if (system == 'Windows') {
        stopCluster(cl)
      }

      return(wObjects)
    }

    # any other blocks - executed any number of times between first and final block
    if (l>1 && l<blocks) {
      message(paste('Processing ', blocksize, ' formulae, batch ',
        l, ' of ', blocks, sep=''))
      message(paste('-- index1: ', (1+(blocksize*(l-1))), '; ',
        'index2: ', (blocksize*l), sep=''))

      df <- data[,c(
        which(colnames(data) %in% terms),
        startIndex + ((1+(blocksize*(l-1))):(blocksize*l)))]

      if (nestedParallel == TRUE) {
        if (system == 'Windows') {
          models <- parLapply(cluster,
            formula.list[(1+(blocksize*(l-1))):(blocksize*l)],
            function(f) summary(FUN(formula = f, data = df))$coefficients)
        } else {
          models <- mclapply(
            formula.list[(1+(blocksize*(l-1))):(blocksize*l)],
            function(f) summary(FUN(formula = f, data = df))$coefficients)
        } 
      } else if (nestedParallel == FALSE) {
        models <- lapply(
          formula.list[(1+(blocksize*(l-1))):(blocksize*l)],
          function(f) summary(FUN(formula = f, data = df))$coefficients)
      } else {
        stop('Invalid value for argument nestedParallel. Must be TRUE/FALSE')
      }

      names(models) <- variables[(1+(blocksize*(l-1))):(blocksize*l)]

      # detect failed models (return 'try' error)
      nullindices <- c(
        which(mapply(inherits, models, 'try-error'))
      )
      if (length(nullindices) > 1) {
        print(models[[nullindices[1]]][1])
        stop(paste('\nOne or more models have failed. ',
          'Please recheck your data and model formula. ',
          'In particular, check that your variable names ',
          'are valid, i.e., no spaces, do not begin with ',
          'a number, no hyphens or other punctuation, et ',
          'cetera.', sep=''))
      }

      # convert to data frames
      if (system == 'Windows') {
        wObjects <- parLapply(cluster, names(models),
          function(x) data.frame(
            rep(x, length(rownames(models[[x]]))),
            rownames(models[[x]]),
            models[[x]],
            row.names=rownames(models[[x]])))
      } else {
        wObjects <- mclapply(names(models),
          function(x) data.frame(
            rep(x, length(rownames(models[[x]]))),
            rownames(models[[x]]),
            models[[x]],
            row.names=rownames(models[[x]])))
      }

      # remove intercept and / or specified terms from the output
      if (!is.null(excludeTerms) && excludeIntercept == TRUE) {
        if (system == 'Windows') {
          wObjects <- parLapply(cluster, wObjects,
            function(x) x[grep(
              paste(c('\\(Intercept\\)', excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        } else {
          wObjects <- mclapply(wObjects,
            function(x) x[grep(
              paste(c('\\(Intercept\\)', excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        }
      } else if (is.null(excludeTerms) && excludeIntercept == TRUE) {
        if (system == 'Windows') {
          wObjects <- parLapply(cluster, wObjects,
            function(x) x[grep(
              '\\(Intercept\\)',
              rownames(x),
              invert=TRUE),])
        } else {
          wObjects <- mclapply(wObjects,
            function(x) x[grep(
              '\\(Intercept\\)',
              rownames(x),
              invert=TRUE),])
        }
      } else if (!is.null(excludeTerms) && excludeIntercept == FALSE) {
        if (system == 'Windows') {
          wObjects <- parLapply(cluster, wObjects,
            function(x) x[grep(
              paste(c(excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        } else {
          wObjects <- mclapply(wObjects,
            function(x) x[grep(
              paste(c(excludeTerms), collapse='|'),
              rownames(x),
              invert=TRUE),])
        }
      }

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), stringsAsFactors=FALSE)

      # set colnames
      # detect if the y variable in a glm is a continous variable
      # (model will contain t stat; otherwise, Z)
      if (any(grepl('^t\\.value$', colnames(wObjects))) == TRUE) {
        colnames(wObjects) <- c('Variable', 'Term', 'Beta',
          'StandardError', 't', 'P')
      } else {
        colnames(wObjects) <- c('Variable', 'Term', 'Beta',
          'StandardError', 'Z', 'P')
      }

      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(
        as.character(wObjects$StandardError))
      if (any(grepl('^t$', colnames(wObjects))) == TRUE) {
        wObjects$t <- as.numeric(as.character(wObjects$t))
      } else {
        wObjects$Z <- as.numeric(as.character(wObjects$Z))
      }
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate OR and confidence intervals
      wObjects$OR <- exp(wObjects[,'Beta'])
      wObjects$ORlower <- wObjects$OR * exp(
        qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$ORupper <- wObjects$OR * exp(
        qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      # change Inf and NaN to NA
      if (any(grepl('^t$', colnames(wObjects)) == TRUE)) {
        wObjects$t[is.infinite(wObjects$t)] <- NA
        wObjects$t[wObjects$t == 'NaN'] <- NA
      } else {
        wObjects$Z[is.infinite(wObjects$Z)] <- NA
        wObjects$Z[wObjects$Z == 'NaN'] <- NA
      }
      wObjects$P[is.infinite(wObjects$P)] <- NA
      wObjects$P[wObjects$P == 'NaN'] <- NA
      wObjects$OR[is.infinite(wObjects$OR)] <- NA
      wObjects$OR[wObjects$OR == 'NaN'] <- NA
      wObjects$ORlower[is.infinite(wObjects$ORlower)] <- NA
      wObjects$ORlower[wObjects$ORlower == 'NaN'] <- NA
      wObjects$ORupper[is.infinite(wObjects$ORupper)] <- NA
      wObjects$ORupper[wObjects$ORupper == 'NaN'] <- NA

      return(wObjects)
    }
  }
}
