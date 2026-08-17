
# DATA TO USE
g <- test_data$vis
EX <- getExpression(g, values = "normalized", output = "matrix")
CL <- getCellMetadata(g, output = "data.table")$leiden_clus
LV <- levels(droplevels(factor(CL)))

.fs <- function(...) {
    analyzeData(EX, analyzeParam("feat_stats"), ...)
}

test_that("grouped feat_stats reproduces the aggregate helpers exactly", {
    # Independent ground truth: the GiottoClass helpers this replaces.
    st <- .fs(groups = CL, stats = c("sum", "nnz"))

    avg <- create_average_DT(g,
        meta_data_name = "leiden_clus", expression_values = "normalized"
    )
    det <- create_average_detection_DT(g,
        meta_data_name = "leiden_clus", expression_values = "normalized",
        detection_threshold = 0
    )
    A <- as.matrix(avg)[, paste0("cluster_", LV), drop = FALSE]
    D <- as.matrix(det)[, paste0("cluster_", LV), drop = FALSE]

    # bit-identical, not merely equal: the gini path downstream has rank ties
    # and strict thresholds that a last-bit difference could flip
    expect_identical(st$mean_expr, as.numeric(A))
    expect_equal(st$perc_cells / 100, as.numeric(D), ignore_attr = TRUE)
})

test_that("grouped feat_stats emits the full feats x groups cross product", {
    st <- .fs(groups = CL, stats = "sum")
    expect_equal(nrow(st), nrow(EX) * length(LV))
    # a feature absent from a group is a zero row, not a missing one -- gini is
    # taken over the length-G vector and a gap would change the coefficient
    expect_equal(st[, .N, by = feats][, unique(N)], length(LV))
    expect_false(any(is.na(st$mean_expr)))
})

test_that("groups vary slowest, feats cycling within", {
    # matches the streaming backend's emission order, so a matrix() reshape
    # lands without transposing
    st <- .fs(groups = CL, stats = "sum")
    expect_equal(st$feats[seq_len(nrow(EX))], rownames(EX))
    expect_equal(unique(st$group[seq_len(nrow(EX))]), LV[1])
})

test_that("accumulator selection drives which columns appear", {
    expect_named(.fs(groups = CL, stats = "sum"),
        c("feats", "group", "n_cells", "total_expr", "mean_expr"))
    expect_named(.fs(groups = CL, stats = "nnz"),
        c("feats", "group", "n_cells", "nr_cells", "perc_cells"))
    expect_true(all(
        c("sumsq", "sd") %in% names(.fs(groups = CL, stats = c("sum", "sumsq")))
    ))
    expect_true(
        "mean_expr_det" %in% names(.fs(groups = CL, stats = c("nnz", "sum_det")))
    )
})

test_that("NA groups are excluded from the population", {
    grp <- CL
    grp[seq_len(50)] <- NA
    st <- .fs(groups = grp, stats = "sum")
    expect_equal(sum(unique(st[, .(group, n_cells)])$n_cells), length(CL) - 50L)
})

test_that("the detection threshold gates counts but never the sums", {
    # adr/0009 semantics, matching the streaming backend. The threshold must
    # sit above the smallest stored value to bite at all -- normalized visium
    # has no nonzero entry below ~1.07, so a threshold of 1 excludes nothing.
    thr <- stats::median(EX@x)
    lo <- .fs(groups = CL, stats = c("sum", "nnz"))
    hi <- analyzeData(EX, analyzeParam("feat_stats", detection_threshold = thr),
        groups = CL, stats = c("sum", "nnz"))
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
    expect_error(.fs(groups = CL[-1L]), "one entry per column")
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
    expect_false(isTRUE(all.equal(
        got$mean_expr, .fs(groups = unname(CL), stats = "sum")$mean_expr
    )))
})

test_that("an unnamed groups vector is still taken positionally", {
    # GiottoDisk passes an unnamed vector aligned to the current view; that
    # contract has to keep working.
    expect_silent(.fs(groups = unname(CL), stats = "sum"))
})

test_that("a named groups vector missing cells errors", {
    ids <- getCellMetadata(g, output = "data.table")[["cell_ID"]]
    named <- stats::setNames(CL, ids)
    expect_error(.fs(groups = named[-1L]), "does not cover every column")
})
