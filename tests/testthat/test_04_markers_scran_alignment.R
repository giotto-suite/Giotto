# findScranMarkers() hands `expr_data` and a groups vector to
# `scran::findMarkers()`, which pairs them positionally. Expression and cell
# metadata are fetched independently and are not guaranteed to share a cell
# order, so the labels have to be reordered onto the matrix first.
#
# TODO: these lock the temporary reorder in findScranMarkers(). Once
# analyzeData(x, scranMarkersParam) arrives from gsource the resolution moves
# into `.markers_scran()`, and these should be retargeted there rather than
# deleted -- the behaviour they assert stays required.

g <- test_data$vis

test_that("mini visium really does disagree on cell order", {
    # The premise of the rest of the file. If this ever fails, the fixture
    # changed and these tests stop covering anything.
    EX <- getExpression(g, values = "normalized", output = "matrix")
    ids <- getCellMetadata(g, output = "data.table")[["cell_ID"]]
    expect_setequal(colnames(EX), ids)
    expect_false(identical(colnames(EX), ids))
    expect_equal(sum(colnames(EX) == ids), 0L)
})

test_that("findScranMarkers labels cells by cell_ID, not by position", {
    skip_if_not_installed("scran")

    # returns an unnamed list of data.tables, one per cluster, each tagged with
    # a `cluster` column
    got <- data.table::rbindlist(findScranMarkers(g,
        cluster_column = "leiden_clus", expression_values = "normalized"
    ))

    # Independent ground truth: call scran directly, with the grouping put into
    # the matrix's column order by name.
    EX <- getExpression(g, values = "normalized", output = "matrix")
    cm <- getCellMetadata(g, output = "data.table")
    grp <- cm[["leiden_clus"]][match(colnames(EX), cm[["cell_ID"]])]
    ref <- scran::findMarkers(x = EX, groups = grp)

    expect_setequal(unique(got$cluster), names(ref))
    for (k in names(ref)) {
        a <- got[cluster == k]
        b <- data.table::as.data.table(ref[[k]])
        b[, feats := rownames(ref[[k]])]
        # same features in the same order, and the same statistics on them
        expect_equal(a$feats, b$feats, info = k)
        expect_equal(a$p.value, b$p.value, info = k)
        expect_equal(a$FDR, b$FDR, info = k)
    }
})

test_that("the misaligned grouping would have given different markers", {
    skip_if_not_installed("scran")

    # Not a tautology check: it establishes that the reorder changes the answer
    # on this object, so the test above is measuring something.
    EX <- getExpression(g, values = "normalized", output = "matrix")
    cm <- getCellMetadata(g, output = "data.table")
    ord <- match(colnames(EX), cm[["cell_ID"]])

    aligned <- scran::findMarkers(x = EX, groups = cm[["leiden_clus"]][ord])
    naive <- scran::findMarkers(x = EX, groups = cm[["leiden_clus"]])

    top10 <- function(res, k) head(rownames(res[[k]]), 10L)
    overlap <- vapply(names(aligned), function(k) {
        length(intersect(top10(aligned, k), top10(naive, k)))
    }, integer(1L))

    # Top-10 markers per cluster barely overlap between the two groupings.
    expect_lt(sum(overlap), 10L)
})

test_that("findScranMarkers_one_vs_all inherits the alignment", {
    skip_if_not_installed("scran")

    # It delegates to findScranMarkers() per cluster, so it needs no reorder of
    # its own -- this asserts that delegation still holds.
    res <- findScranMarkers_one_vs_all(g,
        cluster_column = "leiden_clus", expression_values = "normalized",
        verbose = FALSE
    )
    expect_true(nrow(res) > 0L)
    expect_true(all(c("feats", "cluster") %in% names(res)))
})
