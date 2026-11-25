test_that("library normalization scales to the median library size", {
    m <- matrix(c(1, 3, 0, 4, 2, 2, 0, 0, 5, 0, 0, 0), nrow = 3,
        dimnames = list(paste0("f", 1:3), paste0("c", 1:4)))
    lib <- colSums(m)  # 4 8 5 0; the empty cell is left out of the median
    keep <- lib > 0
    expected <- sweep(m[, keep], 2, lib[keep], "/") * 5
    p <- normParam("library", scalefactor = "median")
    res <- suppressWarnings(processData(m, p))
    expect_equal(res[, keep], expected)
    resd <- suppressWarnings(processData(Matrix::Matrix(m, sparse = TRUE), p))
    expect_equal(as.matrix(resd)[, keep], expected)
    # numeric scalefactors are unchanged; other strings are rejected
    expect_equal(processData(m[, keep], normParam("library",
        scalefactor = 100)), sweep(m[, keep], 2, lib[keep], "/") * 100)
    expect_error(processData(m[, keep], normParam("library",
        scalefactor = "mean")))
})
