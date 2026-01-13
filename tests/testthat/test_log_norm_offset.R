test_that(".log_norm_giotto rejects offset != 1 for sparse matrices", {
    skip_if_not_installed("Matrix")

    mat <- Matrix::rsparsematrix(
        nrow = 10,
        ncol = 8,
        density = 0.15,
        rand.x = function(n) sample.int(10L, n, replace = TRUE) - 1L
    )

    expect_error(
        Giotto:::.log_norm_giotto(mymatrix = mat, base = 2, offset = 0.5),
        regexp = "offset != 1"
    )

    out <- Giotto:::.log_norm_giotto(mymatrix = mat, base = 2, offset = 1)
    expect_equal(out@x, log1p(mat@x) / log(2))
})

test_that(".log_norm_giotto produces equivalent results across sparse types", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("DelayedArray")

    mat <- Matrix::rsparsematrix(
        nrow = 10,
        ncol = 8,
        density = 0.15,
        rand.x = function(n) sample.int(10L, n, replace = TRUE) - 1L
    )
    da_sparse <- DelayedArray::DelayedArray(mat)

    log_norm <- Giotto:::.log_norm_giotto
    res_mat <- log_norm(mat, base = 2, offset = 1)
    res_da <- log_norm(da_sparse, base = 2, offset = 1)

    expect_equal(as.matrix(res_mat), as.matrix(res_da), ignore_attr = TRUE)
})

test_that(".log_norm_giotto rejects offset != 1 for sparse DelayedArray", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("DelayedArray")

    mat <- Matrix::rsparsematrix(
        nrow = 10,
        ncol = 8,
        density = 0.15,
        rand.x = function(n) sample.int(10L, n, replace = TRUE) - 1L
    )
    da_sparse <- DelayedArray::DelayedArray(mat)

    expect_error(
        Giotto:::.log_norm_giotto(mymatrix = da_sparse, base = 2, offset = 0.5),
        regexp = "offset != 1"
    )
})

test_that("processData logNormParam produces equivalent results", {
    skip_if_not_installed("Matrix")
    skip_if_not_installed("DelayedArray")

    mat <- Matrix::rsparsematrix(
        nrow = 10,
        ncol = 8,
        density = 0.15,
        rand.x = function(n) sample.int(10L, n, replace = TRUE) - 1L
    )
    da_sparse <- DelayedArray::DelayedArray(mat)

    param <- Giotto::normParam("log", base = 2, offset = 1)

    res_mat <- GiottoClass::processData(mat, param)
    res_da <- GiottoClass::processData(da_sparse, param)

    expect_equal(as.matrix(res_mat), as.matrix(res_da), ignore_attr = TRUE)
})

test_that("processData logNormParam rejects offset != 1 for sparse", {
    skip_if_not_installed("Matrix")

    mat <- Matrix::rsparsematrix(
        nrow = 10,
        ncol = 8,
        density = 0.15,
        rand.x = function(n) sample.int(10L, n, replace = TRUE) - 1L
    )

    param <- Giotto::normParam("log", base = 2, offset = 0.5)

    expect_error(GiottoClass::processData(mat, param), regexp = "offset != 1")
})

test_that("processData logNormParam allows offset != 1 for dense Matrix", {
    skip_if_not_installed("Matrix")

    # dgeMatrix is a dense Matrix - offset should be allowed
    mat_dense <- Matrix::Matrix(
        matrix(rpois(80, 5), nrow = 10, ncol = 8),
        sparse = FALSE
    )
    expect_true(inherits(mat_dense, "denseMatrix"))

    param <- Giotto::normParam("log", base = 2, offset = 0.5)

    # Should NOT error - dense matrices support any offset
    res <- GiottoClass::processData(mat_dense, param)
    expect_true(inherits(res, "Matrix"))

    # Verify calculation is correct
    expected <- log(as.matrix(mat_dense) + 0.5) / log(2)
    expect_equal(as.matrix(res), expected, ignore_attr = TRUE)
})
