
# DATA TO USE
g <- test_data$vis
EX <- getExpression(g, values = "normalized", output = "matrix")
IDS <- getCellMetadata(g, output = "data.table")[["cell_ID"]]
CL <- getCellMetadata(g, output = "data.table")$leiden_clus
LV <- levels(droplevels(factor(CL)))

# Named by cell ID, which is how every caller should pass a grouping and the
# only form that is correct on this object -- its metadata and expression
# disagree on cell order. `CL` stays unnamed for the tests about that path.
CLN <- stats::setNames(CL, IDS)

.fs <- function(...) {
    analyzeData(EX, analyzeParam("feat_stats"), ...)
}

# Cells of one group, resolved by ID rather than position.
.cells_in <- function(k, groups = CLN) colnames(EX) %in% IDS[groups == k]

test_that("grouped feat_stats reproduces the per-group aggregates exactly", {
    # Ground truth computed here rather than taken from GiottoClass'
    # create_average_DT: that helper was the original home of this statistic
    # and shared the alignment bug, so it cannot referee the fix. It is also a
    # moving target across GiottoClass versions.
    st <- .fs(groups = CLN, stats = c("sum", "nnz"))

    A <- as.numeric(vapply(LV, function(k) {
        Matrix::rowMeans(EX[, .cells_in(k), drop = FALSE])
    }, numeric(nrow(EX))))
    D <- as.numeric(vapply(LV, function(k) {
        Matrix::rowMeans(EX[, .cells_in(k), drop = FALSE] > 0)
    }, numeric(nrow(EX))))

    # bit-identical, not merely equal: the gini path downstream has rank ties
    # and strict thresholds that a last-bit difference could flip
    expect_identical(st$mean_expr, A)
    expect_equal(st$perc_cells / 100, D, ignore_attr = TRUE)
})

test_that("grouped feat_stats emits the full feats x groups cross product", {
    st <- .fs(groups = CLN, stats = "sum")
    expect_equal(nrow(st), nrow(EX) * length(LV))
    # a feature absent from a group is a zero row, not a missing one -- gini is
    # taken over the length-G vector and a gap would change the coefficient
    expect_equal(st[, .N, by = feats][, unique(N)], length(LV))
    expect_false(any(is.na(st$mean_expr)))
})

test_that("groups vary slowest, feats cycling within", {
    # matches the streaming backend's emission order, so a matrix() reshape
    # lands without transposing
    st <- .fs(groups = CLN, stats = "sum")
    expect_equal(st$feats[seq_len(nrow(EX))], rownames(EX))
    expect_equal(unique(st$group[seq_len(nrow(EX))]), LV[1])
})

test_that("accumulator selection drives which columns appear", {
    expect_named(.fs(groups = CLN, stats = "sum"),
        c("feats", "group", "n_cells", "total_expr", "mean_expr"))
    expect_named(.fs(groups = CLN, stats = "nnz"),
        c("feats", "group", "n_cells", "nr_cells", "perc_cells"))
    expect_true(all(
        c("sumsq", "sd") %in% names(.fs(groups = CLN, stats = c("sum", "sumsq")))
    ))
    expect_true(
        "mean_expr_det" %in% names(.fs(groups = CLN, stats = c("nnz", "sum_det")))
    )
})

test_that("NA groups are excluded from the population", {
    grp <- CLN
    grp[seq_len(50)] <- NA
    st <- .fs(groups = grp, stats = "sum")
    expect_equal(sum(unique(st[, .(group, n_cells)])$n_cells), length(CL) - 50L)
})

test_that("the detection threshold gates counts but never the sums", {
    # adr/0009 semantics, matching the streaming backend. The threshold must
    # sit above the smallest stored value to bite at all -- normalized visium
    # has no nonzero entry below ~1.07, so a threshold of 1 excludes nothing.
    thr <- stats::median(EX@x)
    lo <- .fs(groups = CLN, stats = c("sum", "nnz"))
    hi <- analyzeData(EX, analyzeParam("feat_stats", detection_threshold = thr),
        groups = CLN, stats = c("sum", "nnz"))
    expect_equal(lo$total_expr, hi$total_expr)
    expect_lt(sum(hi$nr_cells), sum(lo$nr_cells))
})

test_that("the ungrouped path is unchanged and refuses stats selection", {
    expect_named(.fs(), c(
        "feats", "nr_cells", "perc_cells", "total_expr",
        "mean_expr", "mean_expr_det"
    ))
    expect_error(.fs(stats = "sum"), "grouped path")
})

test_that("an unnamed groups vector of the wrong length errors", {
    # unnamed, so it must be the length check that fires and not the ID match
    expect_error(.fs(groups = unname(CL)[-1L]), "one entry per column")
})


# alignment ####

test_that("a named groups vector is matched on cell ID, not position", {
    # This object's metadata and expression disagree on cell order, which is
    # exactly the case a positional vector gets silently wrong.
    ids <- getCellMetadata(g, output = "data.table")[["cell_ID"]]
    expect_false(identical(ids, colnames(EX)))

    named <- stats::setNames(CL, ids)
    got <- .fs(groups = named, stats = "sum")

    ref <- as.numeric(vapply(LV, function(k) {
        Matrix::rowMeans(EX[, colnames(EX) %in% ids[CL == k], drop = FALSE])
    }, numeric(nrow(EX))))
    expect_equal(got$mean_expr, ref)

    # and the naming genuinely matters here
    pos <- suppressWarnings(.fs(groups = unname(CL), stats = "sum"))
    expect_false(isTRUE(all.equal(got$mean_expr, pos$mean_expr)))
})

test_that("an unnamed groups vector is still taken positionally, with a warning", {
    # A caller that really does hold a column-ordered vector keeps working; the
    # warning is there because nothing can check that claim.
    expect_warning(.fs(groups = unname(CL), stats = "sum"), "unnamed")
})

test_that("cells a named groups vector omits drop out of the statistics", {
    ids <- getCellMetadata(g, output = "data.table")[["cell_ID"]]
    named <- stats::setNames(CL, ids)

    # Dropping every cell of one cluster should remove that group entirely,
    # rather than error -- masking is how callers narrow.
    drop <- names(named)[named == LV[[1L]]]
    got <- .fs(groups = named[!names(named) %in% drop], stats = "sum")
    expect_false(LV[[1L]] %in% got$group)
    expect_setequal(got$group, setdiff(LV, LV[[1L]]))

    # A surviving group's statistics are untouched by the omission.
    ref <- .fs(groups = named, stats = "sum")
    keep <- LV[[2L]]
    expect_equal(
        got[got$group == keep, ]$mean_expr, ref[ref$group == keep, ]$mean_expr
    )

    # Names that match nothing at all is a mistake, not an empty selection.
    expect_error(
        .fs(groups = stats::setNames(CL, paste0("nope_", ids))),
        "none of its names"
    )
})
