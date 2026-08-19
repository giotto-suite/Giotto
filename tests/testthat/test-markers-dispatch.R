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


test_that("non-t tests do not receive the t-only std.lfc argument", {
    skip_if_not_installed("scran")
    fx <- .mk_gobject()

    # `markersParam()` always materializes `std_lfc`, so forwarding it blindly
    # broke every non-t test with "unused argument (std.lfc = FALSE)" --
    # pairwiseWilcox()/pairwiseBinom() have no such parameter. This is the
    # in-memory fallback the streaming backend points users at, so it has to
    # work.
    for (tt in c("wilcox", "binom")) {
        res <- analyzeData(as.matrix(fx$mat),
            markersParam(method = "scran", test_type = tt),
            groups = fx$clus)
        expect_named(as.list(res), c("a", "b", "c"), info = tt)
    }

    # ... and asking for it explicitly is an error rather than a silent drop:
    # an AUC has no standardized-effect-size counterpart.
    expect_error(
        analyzeData(as.matrix(fx$mat),
            markersParam(method = "scran", test_type = "wilcox",
                std_lfc = TRUE),
            groups = fx$clus),
        "applies only to"
    )
})


# `cluster` column type invariant across marker methods. findMarkers_one_vs_all()
# dispatches to all three, so a caller joining or rbind-ing results must be able
# to rely on one type regardless of `method`. A numeric cluster column is used
# deliberately -- that is where a missing coercion shows up, and cluster columns
# from Leiden are numeric.

.mk_numeric_clus_gobject <- function(n_genes = 60L, n_cells = 120L, seed = 1L) {
    set.seed(seed)
    cells <- paste0("c_", seq_len(n_cells))
    m <- matrix(rpois(n_genes * n_cells, 3L), nrow = n_genes,
        dimnames = list(paste0("g_", seq_len(n_genes)), cells))
    clus <- rep(1:3, length.out = n_cells)
    # separate the groups so every method finds markers
    for (k in 1:3) {
        rows <- ((k - 1L) * 10L + 1L):(k * 10L)
        m[rows, clus == k] <- m[rows, clus == k] + 25L
    }
    g <- GiottoClass::createGiottoObject(expression = m, verbose = FALSE)
    g <- Giotto::normalizeGiotto(g, verbose = FALSE)
    GiottoClass::addCellMetadata(g,
        new_metadata = data.frame(cell_ID = cells, clus = clus))
}

test_that("cluster column is character for every marker method", {
    skip_if_not_installed("scran")
    g <- .mk_numeric_clus_gobject()
    # the fixture must actually be numeric, or this proves nothing
    expect_true(is.numeric(GiottoClass::pDataDT(g)$clus))

    scran_ova <- findScranMarkers_one_vs_all(g, cluster_column = "clus",
        verbose = FALSE)
    gini_ova <- findGiniMarkers_one_vs_all(g, cluster_column = "clus",
        verbose = FALSE)

    expect_type(scran_ova$cluster, "character")
    expect_type(gini_ova$cluster, "character")

    # plain (non one-vs-all) paths emit `cluster` too
    scran_plain <- data.table::rbindlist(
        findScranMarkers(g, cluster_column = "clus", verbose = FALSE),
        fill = TRUE)
    expect_type(scran_plain$cluster, "character")
    expect_type(findGiniMarkers(g, cluster_column = "clus")$cluster, "character")

    # the dispatcher must not vary the type by method
    for (m in c("scran", "gini")) {
        res <- findMarkers_one_vs_all(g, method = m, cluster_column = "clus",
            verbose = FALSE)
        expect_type(res$cluster, "character")
    }

    # and the results must be joinable on it -- the reason the type matters
    j <- merge(gini_ova[, c("feats", "cluster")],
               scran_ova[, c("feats", "cluster")],
               by = c("feats", "cluster"))
    expect_gt(nrow(j), 0L)
})

test_that("MAST labels comparisons rather than bare clusters", {
    skip_if_not_installed("MAST")
    g <- .mk_numeric_clus_gobject()
    res <- findMastMarkers_one_vs_all(g, cluster_column = "clus",
        verbose = FALSE)
    # character like the others, but the value is "<clus>_vs_others", not
    # "<clus>" -- so cross-method joins on `cluster` need a translation step.
    expect_type(res$cluster, "character")
    expect_true(all(grepl("_vs_others$", res$cluster)))
})
