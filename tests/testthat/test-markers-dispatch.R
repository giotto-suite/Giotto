# findScranMarkers now routes through analyzeData dispatch rather than calling
# scran directly. The point of these tests is that in memory this changed
# nothing: the result must stay bit-identical to calling scran::findMarkers on
# the same matrix. If a future edit quietly moves the in-memory path onto a
# reimplementation, this is what catches it.

.mk_gobject <- function(n_genes = 20L, n_cells = 45L, seed = 3L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.6,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    clus <- rep(c("a", "b", "c"), length.out = n_cells)
    g <- GiottoClass::createGiottoObject(expression = m)
    g <- GiottoClass::addCellMetadata(g,
        new_metadata = data.frame(clus = clus))
    list(gobject = g, mat = m, clus = clus)
}


test_that("findScranMarkers in memory is identical to scran::findMarkers", {
    skip_if_not_installed("scran")
    fx <- .mk_gobject()

    got <- findScranMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", verbose = FALSE)
    ref <- scran::findMarkers(as.matrix(fx$mat), groups = fx$clus)

    expect_length(got, length(ref))
    for (dt in got) {
        k <- unique(dt$cluster)
        r <- data.table::as.data.table(ref[[k]])
        expect_equal(dt$p.value, r$p.value)
        expect_equal(dt$FDR, r$FDR)
        expect_equal(dt$Top, r$Top)
    }
})


test_that("findScranMarkers forwards scran arguments through scranMarkersParam", {
    skip_if_not_installed("scran")
    fx <- .mk_gobject()

    got <- findScranMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", verbose = FALSE, pval_type = "all")
    ref <- scran::findMarkers(as.matrix(fx$mat), groups = fx$clus,
        pval.type = "all")

    # pval.type = "all" drops the Top column, so this also confirms the
    # snake_case -> dotted translation actually reached scran.
    expect_false("Top" %in% colnames(got[[1L]]))
    for (dt in got) {
        r <- data.table::as.data.table(ref[[unique(dt$cluster)]])
        expect_equal(dt$p.value, r$p.value)
    }
})


test_that("analyzeData(<matrix>, scranMarkersParam) one_vs_rest pools the rest", {
    skip_if_not_installed("scran")
    fx <- .mk_gobject()
    lv <- c("a", "b", "c")

    got <- analyzeData(as.matrix(fx$mat), markersParam(
        method = "scran", comparison = "one_vs_rest"), groups = fx$clus)
    expect_named(as.list(got), lv)

    for (k in lv) {
        rest <- setdiff(lv, k)
        pooled <- ifelse(fx$clus == k, k, paste0(rest, collapse = "_"))
        ref <- scran::findMarkers(as.matrix(fx$mat), groups = pooled)[[k]]
        expect_equal(got[[k]]$p.value, ref$p.value)
    }
})


test_that("findScranMarkers_one_vs_all returns one block per cluster", {
    skip_if_not_installed("scran")
    fx <- .mk_gobject()

    res <- findScranMarkers_one_vs_all(fx$gobject, cluster_column = "clus",
        expression_values = "raw", verbose = FALSE)

    expect_s3_class(res, "data.table")
    expect_setequal(unique(res$cluster), c("a", "b", "c"))
    expect_true(all(c("logFC", "feats", "p.value", "ranking") %in%
        colnames(res)))
})


test_that("findScranMarkers_one_vs_all honours its logFC threshold", {
    skip_if_not_installed("scran")
    fx <- .mk_gobject()

    # The filter is `(p.value <= pval & logFC >= logFC_thresh) |
    # (ranking <= min_feats)`, so `min_feats` is set low enough that the
    # rank clause cannot rescue a feature and the threshold is what decides.
    strict <- findScranMarkers_one_vs_all(fx$gobject, cluster_column = "clus",
        expression_values = "raw", verbose = FALSE,
        logFC = 0.5, min_feats = 1)
    # `-Inf` is the honest "no threshold" control: `logFC = 0` would still
    # exclude down-regulated features, which is a real filter, not the
    # absence of one -- and the shadowed comparison it replaced admitted
    # everything regardless of sign.
    loose <- findScranMarkers_one_vs_all(fx$gobject, cluster_column = "clus",
        expression_values = "raw", verbose = FALSE,
        logFC = -Inf, min_feats = 1)

    # Nothing below the threshold survives on the p-value clause.
    expect_equal(sum(strict$logFC < 0.5 & strict$ranking > 1), 0L)
    # ... and removing the threshold genuinely admits more, which is what
    # fails if the argument is ever shadowed by the column again.
    expect_gt(nrow(loose), nrow(strict))
})
