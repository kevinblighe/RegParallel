test_that('datatypes', {
    expect_type(RegParallel(data),
        c('S4', 'list'))
    expect_type(RegParallel(formula),
        c('language'))
    expect_type(RegParallel(FUN),
        c('closure'))
    expect_type(RegParallel(FUNtype, variables, p.adjust, excludeTerms),
        c('language', 'character'))
    expect_type(RegParallel(blocksize, cores, conflevel),
        c('integer'))
    expect_gt(RegParallel(blocksize, cores, conflevel),
        0)
    expect_type(RegParallel(nestedParallel, excludeIntercept),
        c('logical'))
})
