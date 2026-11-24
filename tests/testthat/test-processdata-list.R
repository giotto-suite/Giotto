test_that("a list of params runs in order on data without its own list method", {
    # A container with methods for single params only
    setClass("procListBox", slots = c(v = "numeric"))
    setMethod("processData", signature("procListBox", "logNormParam"),
        function(x, param, ...) {
            x@v <- log(x@v + param$offset, base = param$base)
            x
        })
    setMethod("processData", signature("procListBox", "libraryNormParam"),
        function(x, param, ...) {
            x@v <- x@v / sum(x@v) * param$scalefactor
            x
        })
    on.exit({
        removeMethod("processData", signature("procListBox", "logNormParam"))
        removeMethod("processData", signature("procListBox", "libraryNormParam"))
        removeClass("procListBox")
    }, add = TRUE)

    x <- new("procListBox", v = c(1, 3))
    lib <- normParam("library", scalefactor = 100)
    lg <- normParam("log", base = 2, offset = 1)
    res <- processData(x, list(lib, lg))
    expect_equal(res@v, log2(c(25, 75) + 1))
    # Order matters
    expect_equal(processData(x, list(lg, lib))@v, log2(c(2, 4)) / sum(log2(c(2, 4))) * 100)
    expect_identical(processData(x, list())@v, x@v)
    expect_error(processData(x, list(lib, 1)), "not recognized")

    # Classes with their own list method keep it
    expect_identical(selectMethod("processData", c("dgCMatrix", "list"))@defined[[1L]], "allMatrix")
})
