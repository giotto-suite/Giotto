# Atera output is read on disk only (GiottoDisk::importAteraDisk()). There is no
# in-memory Atera reader: both entry points stop without a backend, before
# touching the data. Deliberately data-free so they run in CI.

test_that("importAtera() without a backend stops and points to backend =", {
    expect_error(importAtera(tempdir()), "backend")
    expect_error(importAtera(tempdir()), "importXenium")
})

test_that("createGiottoAteraObject() without a backend stops the same way", {
    expect_error(createGiottoAteraObject(tempdir()), "backend")
})

test_that("there is no in-memory Atera reader class", {
    expect_false(methods::isClass("AteraReader"))
})

test_that("the Xenium reader keeps its own platform label", {
    expect_identical(getClass("XeniumReader")@prototype@platform, "Xenium")
    expect_identical(
        formals(Giotto:::.xenium_detect_paths)$platform, "Xenium"
    )
})
