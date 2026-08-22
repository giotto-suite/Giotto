# Atera output currently has the same layout as Xenium output, so `AteraReader`
# subclasses `XeniumReader` and overrides nothing but the platform label used in
# path-detection messages. These assertions are the contract that keeps that
# true -- if someone later copies the Xenium reader instead of extending it, or
# renames the label slot, these fail. Deliberately data-free so they run in CI.

test_that("AteraReader extends XeniumReader", {
    expect_true(extends("AteraReader", "XeniumReader"))
})

test_that("the platform label is Atera, and Xenium's is unchanged", {
    expect_identical(getClass("AteraReader")@prototype@platform, "Atera")
    expect_identical(getClass("XeniumReader")@prototype@platform, "Xenium")
})

test_that("path detection defaults to the Xenium label", {
    # the default must be preserved for every existing caller
    expect_identical(
        formals(Giotto:::.xenium_detect_paths)$platform, "Xenium"
    )
})

test_that("$ and autocomplete are inherited", {
    expect_identical(
        .DollarNames.AteraReader, .DollarNames.XeniumReader
    )
})
