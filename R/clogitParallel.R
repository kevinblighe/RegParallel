clogitParallel <- function(
  data,
  formula.list,
  FUN,
  variables,
  terms,
  startIndex,
  blocksize,
  blocks,
  APPLYFUN,
  conflevel,
  excludeTerms)
{

  ExpBeta <- l <- NULL

  # loop through and process each block of variants
  # with foreach, this loop is parallelised
  foreach(l = seq_len(blocks),
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table', 'doParallel', 'survival')) %dopar% {

    # first block - will be executed just once
    if (l==1) {
      message('Processing ', blocksize,
        ' formulae, batch ', l, ' of ', blocks)
      message('-- index1: 1; ', 'index2: ', (blocksize*l))

      # subset the data to only include variants that will be
      # processsed in this block
      df <- data[,c(
        which(colnames(data) %in% terms),
        startIndex + (seq_len(blocksize*l)))]

      # run the models
      models <- APPLYFUN(formula.list[seq_len(blocksize*l)],
        function(f) summary(FUN(formula = f, data = df)))
      names(models) <- variables[seq_len(blocksize*l)]

      # detect failed models (return 'try' error)
      nullindices <- c(
        which(mapply(inherits, models, 'try-error'))
      )
      if (length(nullindices) > 1) {
        print(models[[nullindices[1]]][1])
        stop('\nOne or more models have failed. ',
          'Please recheck your data and model formula. ',
          'In particular, check that your variable names ',
          'are valid, i.e., no spaces, do not begin with ',
          'a number, no hyphens or other punctuation, et ',
          'cetera.')
      }

      # for coxph, clogit, and glm.nb, we have to extract 2 types of
      # information from the models:
      #   data at level of coefficients
      #   data at level of model
      # Here, extract coefficients
      wObjects <- APPLYFUN(names(models),
        function(x) data.frame(models[[x]]['coefficients']))
      names(wObjects) <- variables[seq_len(blocksize*l)]

      # further processing on coefficient data (form into a dataframe) as
      # per glm, lm, bayesglm
      wObjects <- APPLYFUN(names(wObjects),
        function(x) data.frame(
          rep(x, length(rownames(wObjects[[x]]))),
          rownames(wObjects[[x]]),
          wObjects[[x]],
          row.names=rownames(wObjects[[x]])))

      # extract model-level data
      p <- APPLYFUN(names(models),
        function(x) data.frame(
          models[[x]]$logtest,
          models[[x]]$waldtest,
          models[[x]]$sctest)['pvalue',])

      # remove specified terms from the output
      # Nota bene - intercept not relevant for Cox models
      if (!is.null(excludeTerms) ) {
        wObjects <- APPLYFUN(wObjects,
          function(x) x[grep(
            paste(c(excludeTerms), collapse='|'),
            rownames(x),
            invert=TRUE),])
      }

      # check the final number of terms per variable and replicate the coxph-
      # specific p-value data-frame rows by these.
      # As there is only one of these p-values returned per model we have to
      # duplicate them to fit alongside the respective model terms in the final data-table
      nterms <- nrow(wObjects[[1]])
      p <- APPLYFUN(p, function(x) x[rep('pvalue', nterms),])

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), rbindlist(p), stringsAsFactors=FALSE)

      # set colnames
      colnames(wObjects) <- c('Variable', 'Term', 'Beta', 'ExpBeta',
        'StandardError', 'Z', 'P', 'LRT', 'Wald', 'LogRank')

      # sort out final data to return
      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(
        as.character(wObjects$StandardError))
      wObjects$Z <- as.numeric(as.character(wObjects$Z))
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate HR and confidence intervals
      #wObjects$HR <- exp(wObjects[,'Beta'])
      wObjects$HR <- wObjects$ExpBeta
      wObjects$HRlower <- wObjects$HR * exp(
        qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$HRupper <- wObjects$HR * exp(
        qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      # change Inf and NaN to NA
      wObjects$Z[is.infinite(wObjects$Z)] <- NA
      wObjects$Z[wObjects$Z == 'NaN'] <- NA
      wObjects$P[is.infinite(wObjects$P)] <- NA
      wObjects$P[wObjects$P == 'NaN'] <- NA
      wObjects$HR[is.infinite(wObjects$HR)] <- NA
      wObjects$HR[wObjects$HR == 'NaN'] <- NA
      wObjects$HRlower[is.infinite(wObjects$HRlower)] <- NA
      wObjects$HRlower[wObjects$HRlower == 'NaN'] <- NA
      wObjects$HRupper[is.infinite(wObjects$HRupper)] <- NA
      wObjects$HRupper[wObjects$HRupper == 'NaN'] <- NA
      
      return(wObjects[,ExpBeta:=NULL])
    }

    # final block - will be executed just once
    if (l==blocks) {
      message('Processing final batch ', l, ' of ', blocks)
      message('-- index1: ', (1+(blocksize*(l-1))), '; ',
        'index2: ', length(formula.list))

      # subset the data to only include variants that will be
      # processsed in this block
      df <- data[,c(
        which(colnames(data) %in% terms),
        startIndex + (((1+(blocksize*(l-1)))):(length(formula.list))))]

      # run the models
      models <- APPLYFUN(
        formula.list[(1+(blocksize*(l-1))):length(formula.list)],
        function(f) summary(FUN(formula = f, data = df)))
      names(models) <- variables[(1+(blocksize*(l-1))):length(formula.list)]

      # detect failed models (return 'try' error)
      nullindices <- c(
        which(mapply(inherits, models, 'try-error'))
      )
      if (length(nullindices) > 1) {
        print(models[[nullindices[1]]][1])
        stop('\nOne or more models have failed. ',
          'Please recheck your data and model formula. ',
          'In particular, check that your variable names ',
          'are valid, i.e., no spaces, do not begin with ',
          'a number, no hyphens or other punctuation, et ',
          'cetera.')
      }

      # for coxph, clogit, and glm.nb, we have to extract 2 types of
      # information from the models:
      #   data at level of coefficients
      #   data at level of model
      # Here, extract coefficients
      wObjects <- APPLYFUN(names(models),
        function(x) data.frame(models[[x]]['coefficients']))
      names(wObjects) <- variables[(1+(blocksize*(l-1))):length(formula.list)]

      # further processing on coefficient data (form into a dataframe) as
      # per glm, lm, bayesglm
      wObjects <- APPLYFUN(names(wObjects),
        function(x) data.frame(
          rep(x, length(rownames(wObjects[[x]]))),
          rownames(wObjects[[x]]),
          wObjects[[x]],
          row.names=rownames(wObjects[[x]])))

      # extract model-level data
      p <- APPLYFUN(names(models),
        function(x) data.frame(
          models[[x]]$logtest,
          models[[x]]$waldtest,
          models[[x]]$sctest)['pvalue',])

      # remove specified terms from the output
      # Nota bene - intercept not relevant for Cox models
      if (!is.null(excludeTerms) ) {
        wObjects <- APPLYFUN(wObjects,
          function(x) x[grep(
            paste(c(excludeTerms), collapse='|'),
            rownames(x),
            invert=TRUE),])
      }

      # check the final number of terms per variable and replicate the coxph-
      # specific p-value data-frame rows by these.
      # As there is only one of these p-values returned per model we have to
      # duplicate them to fit alongside the respective model terms in the final data-table
      nterms <- nrow(wObjects[[1]])
      p <- APPLYFUN(p, function(x) x[rep('pvalue', nterms),])

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), rbindlist(p), stringsAsFactors=FALSE)

      # set colnames
      colnames(wObjects) <- c('Variable', 'Term', 'Beta', 'ExpBeta',
        'StandardError', 'Z', 'P', 'LRT', 'Wald', 'LogRank')

      # sort out final data to return
      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(
        as.character(wObjects$StandardError))
      wObjects$Z <- as.numeric(as.character(wObjects$Z))
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate HR and confidence intervals
      #wObjects$HR <- exp(wObjects[,'Beta'])
      wObjects$HR <- wObjects$ExpBeta
      wObjects$HRlower <- wObjects$HR * exp(
        qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$HRupper <- wObjects$HR * exp(
        qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      # change Inf and NaN to NA
      wObjects$Z[is.infinite(wObjects$Z)] <- NA
      wObjects$Z[wObjects$Z == 'NaN'] <- NA
      wObjects$P[is.infinite(wObjects$P)] <- NA
      wObjects$P[wObjects$P == 'NaN'] <- NA
      wObjects$HR[is.infinite(wObjects$HR)] <- NA
      wObjects$HR[wObjects$HR == 'NaN'] <- NA
      wObjects$HRlower[is.infinite(wObjects$HRlower)] <- NA
      wObjects$HRlower[wObjects$HRlower == 'NaN'] <- NA
      wObjects$HRupper[is.infinite(wObjects$HRupper)] <- NA
      wObjects$HRupper[wObjects$HRupper == 'NaN'] <- NA

      return(wObjects[,ExpBeta:=NULL])
    }

    # any other blocks - executed any number of times between first and final block
    if (l>1 && l<blocks) {
      message('Processing ', blocksize, ' formulae, batch ',
        l, ' of ', blocks)
      message('-- index1: ', (1+(blocksize*(l-1))), '; ',
        'index2: ', (blocksize*l))

      # subset the data to only include variants that will be
      # processsed in this block
      df <- data[,c(
        which(colnames(data) %in% terms),
        startIndex + ((1+(blocksize*(l-1))):(blocksize*l)))]

      # run the models
      models <- APPLYFUN(
        formula.list[(1+(blocksize*(l-1))):(blocksize*l)],
        function(f) summary(FUN(formula = f, data = df)))
      names(models) <- variables[(1+(blocksize*(l-1))):(blocksize*l)]

      # detect failed models (return 'try' error)
      nullindices <- c(
        which(mapply(inherits, models, 'try-error'))
      )
      if (length(nullindices) > 1) {
        print(models[[nullindices[1]]][1])
        stop('\nOne or more models have failed. ',
          'Please recheck your data and model formula. ',
          'In particular, check that your variable names ',
          'are valid, i.e., no spaces, do not begin with ',
          'a number, no hyphens or other punctuation, et ',
          'cetera.')
      }

      # for coxph, clogit, and glm.nb, we have to extract 2 types of
      # information from the models:
      #   data at level of coefficients
      #   data at level of model
      # Here, extract coefficients
      wObjects <- APPLYFUN(names(models),
        function(x) data.frame(models[[x]]['coefficients']))
      names(wObjects) <- variables[(1+(blocksize*(l-1))):(blocksize*l)]

      # further processing on coefficient data (form into a dataframe) as
      # per glm, lm, bayesglm
      wObjects <- APPLYFUN(names(wObjects),
        function(x) data.frame(
          rep(x, length(rownames(wObjects[[x]]))),
          rownames(wObjects[[x]]),
          wObjects[[x]],
          row.names=rownames(wObjects[[x]])))

      # extract model-level data
      p <- APPLYFUN(names(models),
        function(x) data.frame(
          models[[x]]$logtest,
          models[[x]]$waldtest,
          models[[x]]$sctest)['pvalue',])

      # remove specified terms from the output
      # Nota bene - intercept not relevant for Cox models
      if (!is.null(excludeTerms) ) {
        wObjects <- APPLYFUN(wObjects,
          function(x) x[grep(
            paste(c(excludeTerms), collapse='|'),
            rownames(x),
            invert=TRUE),])
      }

      # check the final number of terms per variable and replicate the coxph-
      # specific p-value data-frame rows by these.
      # As there is only one of these p-values returned per model we have to
      # duplicate them to fit alongside the respective model terms in the final data-table
      nterms <- nrow(wObjects[[1]])
      p <- APPLYFUN(p, function(x) x[rep('pvalue', nterms),])

      #Convert the list into a data table
      wObjects <- data.table(rbindlist(wObjects), rbindlist(p), stringsAsFactors=FALSE)

      # set colnames
      colnames(wObjects) <- c('Variable', 'Term', 'Beta', 'ExpBeta',
        'StandardError', 'Z', 'P', 'LRT', 'Wald', 'LogRank')

      # sort out final data to return
      wObjects$Variable <- as.character(wObjects$Variable)
      wObjects$Beta <- as.numeric(as.character(wObjects$Beta))
      wObjects$StandardError <- as.numeric(
        as.character(wObjects$StandardError))
      wObjects$Z <- as.numeric(as.character(wObjects$Z))
      wObjects$P <- as.numeric(as.character(wObjects$P))

      # calculate HR and confidence intervals
      #wObjects$HR <- exp(wObjects[,'Beta'])
      wObjects$HR <- wObjects$ExpBeta
      wObjects$HRlower <- wObjects$HR * exp(
        qnorm(((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)
      wObjects$HRupper <- wObjects$HR * exp(
        qnorm(1 - ((1 - (conflevel / 100)) / 2)) * wObjects$StandardError)

      # change Inf and NaN to NA
      wObjects$Z[is.infinite(wObjects$Z)] <- NA
      wObjects$Z[wObjects$Z == 'NaN'] <- NA
      wObjects$P[is.infinite(wObjects$P)] <- NA
      wObjects$P[wObjects$P == 'NaN'] <- NA
      wObjects$HR[is.infinite(wObjects$HR)] <- NA
      wObjects$HR[wObjects$HR == 'NaN'] <- NA
      wObjects$HRlower[is.infinite(wObjects$HRlower)] <- NA
      wObjects$HRlower[wObjects$HRlower == 'NaN'] <- NA
      wObjects$HRupper[is.infinite(wObjects$HRupper)] <- NA
      wObjects$HRupper[wObjects$HRupper == 'NaN'] <- NA

      return(wObjects[,ExpBeta:=NULL])
    }
  }
}
