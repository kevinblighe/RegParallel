test_that("datatypes", {
    expect_type(RegParallel(data),
        c("S4", "list"))
    expect_type(RegParallel(formula),
        c("language"))
    expect_type(RegParallel(variables),
        c("language", "character"))
    expect_type(blocksize, cores),
        c("integer"))
    expect_gt(RegParallel(blocksize, cores),
        0)
})
