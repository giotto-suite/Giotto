test_that(".log_norm_giotto rejects offset != 1 for sparse matrices", {
    skip_if_not_installed("Matrix")

    # Use non-negative values to avoid NaNs from log() on negative numbers
    mat <- Matrix::rsparsematrix(
        nrow = 10,
        ncol = 8,
        density = 0.15,
        rand.x = function(n) sample.int(10L, n, replace = TRUE) - 1L
    )

    expect_error(
        object = Giotto:::.log_norm_giotto(mymatrix = mat, base = 2, offset = 0.5),
        regexp = "offset != 1"
    )

    out <- Giotto:::.log_norm_giotto(mymatrix = mat, base = 2, offset = 1)
    expect_equal(out@x, log(mat@x + 1) / log(2))
})
