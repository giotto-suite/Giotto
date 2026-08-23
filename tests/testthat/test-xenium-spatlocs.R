# Cell centroids as a spatial-locations fallback, used when a Xenium export
# ships no boundary files. The cell metadata file carries x_centroid/y_centroid
# and .xenium_cellmeta() drops them on purpose, so this reads them separately.

.mk_cells_parquet <- function(dir) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    dt <- data.table::data.table(
        cell_id = c("a", "b", "c"),
        x_centroid = c(1.5, 2.5, 3.5),
        y_centroid = c(10.0, 20.0, 30.0),
        total_counts = c(5L, 6L, 7L)
    )
    p <- file.path(dir, "cells.parquet")
    arrow::write_parquet(dt, p)
    p
}

test_that("centroids become spatial locations", {
    skip_if_not_installed("arrow")
    d <- file.path(tempdir(), "xen_sl_1")
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    sl <- Giotto:::.xenium_spatlocs(.mk_cells_parquet(d), verbose = FALSE)

    expect_s4_class(sl, "spatLocsObj")
    dt <- sl[]
    expect_setequal(colnames(dt), c("cell_ID", "sdimx", "sdimy"))
    expect_setequal(dt$cell_ID, c("a", "b", "c"))
    expect_equal(dt[cell_ID == "a"]$sdimx, 1.5)
})

test_that("y is flipped, matching the transcript and polygon loaders", {
    skip_if_not_installed("arrow")
    d <- file.path(tempdir(), "xen_sl_2")
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    p <- .mk_cells_parquet(d)

    # both default to flip_vertical = TRUE, so spatial locations must agree --
    # otherwise cells land mirrored against every other Xenium modality
    expect_true(formals(Giotto:::.xenium_spatlocs)$flip_vertical)
    flipped <- Giotto:::.xenium_spatlocs(p, verbose = FALSE)[]
    expect_equal(flipped[cell_ID == "a"]$sdimy, -10)

    unflipped <- Giotto:::.xenium_spatlocs(p, flip_vertical = FALSE,
        verbose = FALSE)[]
    expect_equal(unflipped[cell_ID == "a"]$sdimy, 10)
})

test_that("an absent cell metadata path yields NULL rather than an error", {
    expect_null(Giotto:::.xenium_spatlocs(character(0)))
    expect_null(Giotto:::.xenium_spatlocs(NULL))
    expect_null(Giotto:::.xenium_spatlocs(""))
})

test_that("the create path falls back to centroids when no polygons load", {
    # guards against the two create_gobject copies drifting: the fallback has
    # to be wired in Giotto's copy as well as GiottoDisk's
    src <- paste(deparse(
        methods::getMethod("initialize", "XeniumReader")@.Data
    ), collapse = " ")
    expect_true(grepl("load_spatlocs", src))
})
