RegParallel <- function(
  data,
  formula,
  FUN,
  FUNtype,
  variables,
  blocksize = 500,
  cores = 3,
  nestedParallel = FALSE,
  conflevel = 95,
  excludeTerms = NULL,
  excludeIntercept = TRUE)
{
  message(paste0('\n##############################\n#',
    'RegParallel\n',
    '##############################\n'))

  # detect system being used
  system <- Sys.info()['sysname']
  message('System is:')
  message('-- ', system)

  # blocksize
  blocksize <- round(blocksize, 0)
  message('Blocksize:')
  message('-- ', blocksize)

  # remove and report data with all 0
  zeros <- which(mapply(function(x) all(x == 0), data[variables]) == TRUE)
  if (length(zeros) > 0 ) {
    message(length(zeros), ' variables have all zeros! Removing these...')
    data <- data[,-which(colnames(data) %in% variables[zeros])]
    variables <- variables[-zeros]
    message('Now testing ', length(variables), ' variables')
  }

  # blocksize cannot be greater than number of variables to test
  # throw an error if so
  if (blocksize > length(variables)) {
    stop(paste('blocksize is greater than number of variables to test.',
      'Choose a smaller blocksize'),
      call. = FALSE)
  }

  # cores / threads
  if (cores < 1) {
    stop('Cannot have cores / threads less than 1')
  }
  message('Cores / Threads:')
  message('-- ', floor(cores))

  # nested parallelisation
  if (nestedParallel == TRUE) {
    message('Nesting enabled. Total potential cores + forked threads is:')
    message('-- ', cores * cores)
  }

  # if user is on Windows, cores / threads for certain functions have
  # to be registered differently, as a 'cl'.
  # If on Mac / Linux, they are registered via 'mc.cores'
  cl <- NULL
  if (system == 'Windows') {
    cl <- makeCluster(getOption('cl.cores', cores))
    registerDoSEQ()
    on.exit(stopCluster(cl))
  } else {
    options('mc.cores' = cores)
  }

  # register cores / threads for doParallel-related functions.
  # works on all platforms
  registerDoParallel(cores)

  # determine number of blocks
  blocks <- floor((length(variables)) / blocksize) + 1
  if (blocksize == length(variables)) {
    blocks <- 1
  }

  #########
  # parsing the formula
    # substitute the [*] for a dummy variable and try to coerce to a formula
    f <- gsub('\\[\\*\\]', 'dummyvar', formula)
    if (is(try(as.formula(f), silent = TRUE), 'try-error')) {
      stop('The formula is invalid.')
    }

    # strip out information from the model to parse the variable names
    terms <- all.vars(as.formula(f))

    # if, on the off chance, the user is testing for a variable called
    # 'dummyvar', we have to account for this possibility
    if (str_count(f, 'dummyvar') > 1) {
      terms <- terms
    } else {
      terms <- terms[-grep('dummyvar', terms)]
    }

    message('Terms included in model:')
    for (i in 1:length(terms)) {
      message('-- ', terms[i])
    }

    # check that terms are in the data
    for (i in 1:length(terms)) {
      if (!any(grepl(paste0('^', terms[i], '$'), colnames(data))) == TRUE) {
        stop(paste('Error! - ', terms[i], ' not found. ',
          'Check your model formula and data.', sep=''))
      }
    }
  #########

  # 'left align' the terms in the input data
  data <- cbind(data[,which(colnames(data) %in% c(terms))], data[,which(colnames(data) %in% c(variables))])

  # determine starting index in the data from which we will loop over our
  # variables being tested
  startIndex <- length(terms)

  # when there is only 1 term, the colname of the data-frame defaults to 
  # 'data[, which(colnames(data) %in% c(terms))]'
  # reset it
  if (length(terms) == 1) {
    colnames(data)[1] <- terms[1]
  }

  # store each possible formula in the list
  formula.list <- list()
  for (i in 1:length(variables)) {
    formula.list[[i]] <- as.formula(gsub('\\[\\*\\]', variables[i], formula))
  }

  # dispaly first 5 formula
  five <- unlist(head(formula.list), 5)
  message("First 5 formulae:")
  for (i in 1:5) {
    message('-- ', five[i])
  }

  # no Intercept term is returned for coxph or clogit
  # issue a message to user regarding this
  if ((excludeIntercept == FALSE) && (grepl('coxph|clogit', FUNtype) == TRUE)) {
    message(paste0('Note: an intercept term will be generated for neither ',
      'Cox Proportional Hazards nor conditional logistic regression ',
      'models.'))
  }

  # detect if incorrect FUNtype selected
  if (!any((FUNtype %in% c('glm', 'lm', 'coxph', 'clogit', 'bayesglm', 'glm.nb')) == TRUE)) {
    stop(paste0('FUNtype not recognised. Choose one of glm, lm, coxph, ',
      'clogit, bayesglm, or glm.nb'))
  }

  # identify FUNtype and launch relevant function
  if (FUNtype == 'glm') {
    res <- glmParallel(
      data = data,
      formula.list = formula.list,
      FUN = FUN,
      variables = variables,
      terms = terms,
      startIndex = startIndex,
      blocksize = blocksize,
      blocks = blocks,
      system = system,
      cl = cl,
      nestedParallel = nestedParallel,
      conflevel = conflevel,
      excludeTerms = excludeTerms,
      excludeIntercept = excludeIntercept)
  } else if (FUNtype == 'lm') {
    res <- lmParallel(
      data = data,
      formula.list = formula.list,
      FUN = FUN,
      variables = variables,
      terms = terms,
      startIndex = startIndex,
      blocksize = blocksize,
      blocks = blocks,
      system = system,
      cl = cl,
      nestedParallel = nestedParallel,
      conflevel = conflevel,
      excludeTerms = excludeTerms,
      excludeIntercept = excludeIntercept)
  } else if (FUNtype == 'coxph') {
    res <- coxphParallel(
      data = data,
      formula.list = formula.list,
      FUN = FUN,
      variables = variables,
      terms = terms,
      startIndex = startIndex,
      blocksize = blocksize,
      blocks = blocks,
      system = system,
      cl = cl,
      nestedParallel = nestedParallel,
      conflevel = conflevel,
      excludeTerms = excludeTerms)
  } else if (FUNtype == 'clogit') {
    res <- clogitParallel(
      data = data,
      formula.list = formula.list,
      FUN = FUN,
      variables = variables,
      terms = terms,
      startIndex = startIndex,
      blocksize = blocksize,
      blocks = blocks,
      system = system,
      cl = cl,
      nestedParallel = nestedParallel,
      conflevel = conflevel,
      excludeTerms = excludeTerms)
  } else if (FUNtype == 'bayesglm') {
    res <- bayesglmParallel(
      data = data,
      formula.list = formula.list,
      FUN = FUN,
      variables = variables,
      terms = terms,
      startIndex = startIndex,
      blocksize = blocksize,
      blocks = blocks,
      system = system,
      cl = cl,
      nestedParallel = nestedParallel,
      conflevel = conflevel,
      excludeTerms = excludeTerms,
      excludeIntercept = excludeIntercept)
  } else if (FUNtype == 'glm.nb') {
    res <- glm.nbParallel(
      data = data,
      formula.list = formula.list,
      FUN = FUN,
      variables = variables,
      terms = terms,
      startIndex = startIndex,
      blocksize = blocksize,
      blocks = blocks,
      system = system,
      cl = cl,
      nestedParallel = nestedParallel,
      conflevel = conflevel,
      excludeTerms = excludeTerms,
      excludeIntercept = excludeIntercept)
  }

  # a sanity / integrity check - this code should never be executed
  # checks that all variables in results object are the same and are in the
  # same order as they input.
  if (!all((unique(res$Variable) == variables) == TRUE)) {
    stop(paste0('An unknown error has occurred. ',
      'Please verify that adequate resources are available and retry'))
  }

  message('Done!')

  # if Windows system, disable access to grabbed cl
  #if (system == 'Windows') {
  #  stopcl(cl)
  #}

  return(res)
}
