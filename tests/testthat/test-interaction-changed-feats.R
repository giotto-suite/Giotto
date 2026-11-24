.icf_fixture <- function(backend = NULL) {
    withr::local_seed(1)
    n <- 300
    cells <- paste0("c", seq_len(n))
    m <- Matrix::Matrix(matrix(stats::rpois(15 * n, 3), 15, n,
        dimnames = list(paste0("g", 1:15), cells)), sparse = TRUE)
    locs <- data.frame(sdimx = stats::runif(n, 0, 100),
        sdimy = stats::runif(n, 0, 100), cell_ID = cells)
    g <- createGiottoObject(expression = m, spatial_locs = locs,
        backend = backend, verbose = FALSE)
    g <- normalizeGiotto(g, scale_feats = FALSE, scale_cells = FALSE,
        verbose = FALSE)
    g <- addCellMetadata(g, new_metadata = data.frame(cell_ID = cells,
        type = sample(c("A", "B", "C"), n, TRUE)),
        by_column = TRUE, column_cell_ID = "cell_ID")
    createSpatialNetwork(g, method = "Delaunay", verbose = FALSE)
}

.icf_scores <- function(g, diff_test = "t.test", ...) {
    res <- findInteractionChangedFeats(g, cluster_column = "type",
        diff_test = diff_test, do_parallel = FALSE, ...)$ICFscores
    data.table::setorder(res, unif_int, cell_type, feats)
    res
}

test_that("selected_feats restricts the features tested", {
    g <- .icf_fixture()
    res <- .icf_scores(g, selected_feats = c("g2", "g5", "absent"))
    expect_setequal(unique(res$feats), c("g2", "g5"))
})

test_that("one selected feature gives one row per interaction", {
    g <- .icf_fixture()
    cols <- c("feats", "unif_int", "cell_type", "int_cell_type", "sel",
        "other", "log2fc", "diff", "p.value")
    for (tst in c("t.test", "wilcox")) {
        one <- .icf_scores(g, diff_test = tst, selected_feats = "g2")
        two <- .icf_scores(g, diff_test = tst, selected_feats = c("g2", "g5"))
        expect_identical(unique(one$feats), "g2")
        expect_equal(as.data.frame(one[, ..cols]),
            as.data.frame(two[feats == "g2", ..cols]))
    }
    ref <- as.data.frame(one[, c("unif_int", "cell_type", "sel", "other")])
    for (tst in c("limma", "permutation")) {
        res <- .icf_scores(g, diff_test = tst, selected_feats = "g2",
            nr_permutations = 20)
        expect_equal(as.data.frame(res[, c("unif_int", "cell_type", "sel",
            "other")]), ref)
        expect_true(all(is.finite(res$p.value)))
    }
})

test_that("a GiottoDisk-backed object gives the in-memory result", {
    skip_if_not_installed("GiottoDisk")
    vault <- withr::local_tempdir()
    mem <- .icf_fixture()
    disk <- .icf_fixture(backend = file.path(vault, "proj"))
    expect_true(inherits(getExpression(disk, values = "normalized")[],
        "dataStore"))
    cols <- c("feats", "cell_type", "int_cell_type", "nr_select", "sel",
        "other", "log2fc", "p.value")
    expect_equal(.icf_scores(disk)[, ..cols], .icf_scores(mem)[, ..cols],
        tolerance = 1e-10)
    expect_equal(.icf_scores(disk, selected_feats = c("g1", "g3"))[, ..cols],
        .icf_scores(mem, selected_feats = c("g1", "g3"))[, ..cols],
        tolerance = 1e-10)
})
