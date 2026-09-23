# Differential expression at every branch point of a cluster tree.
#
# `findNodeMarkers()` has two routes to the same answer: a pooled one, which
# takes a single grouped statistic pass and combines the accumulators per node,
# and a delegated one, which calls `findMarkers(group_1 =, group_2 =)` once per
# node. The pooled route is the fast path and the delegated route is the
# definition, so the load-bearing assertion in this file is that they agree.
#
# The second thing pinned here is the gini parameter defaults. `markersParam()`
# defaults looser than `findMarkers()` does (0.2 / 0.2 and `min_feats` 5, against
# 0.5 / 0.5 and 4). Inheriting those silently made a node comparison a different
# statistic from the pairwise comparison it is meant to match, returning an order
# of magnitude more rows.

.nm_gobject <- function(n_genes = 60L, n_cells = 120L, n_clus = 6L, seed = 11L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.7,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    clus <- as.character(rep(seq_len(n_clus), length.out = n_cells))
    # give each cluster a block of genes it is high in, so the tree has
    # structure rather than being noise
    for (k in seq_len(n_clus)) {
        gi <- ((k - 1L) * 5L + 1L):(k * 5L)
        m[gi, clus == as.character(k)] <- m[gi, clus == as.character(k)] + 20
    }
    g <- GiottoClass::createGiottoObject(expression = m)
    # keyed by cell_ID, not positional: addCellMetadata warns otherwise
    g <- GiottoClass::addCellMetadata(g,
        new_metadata = data.frame(cell_ID = colnames(m), clus = clus))
    list(gobject = g, mat = m, clus = clus)
}

.nm_tree <- function(g) {
    calculateClusterTree(g, cluster_column = "clus",
        expression_values = "raw", cor = "pearson", distance = "average")
}

# the scran result is an UNNAMED list whose elements carry a `cluster` column
# naming their host; pick the host that is this node's left side
.nm_host <- function(res, name) {
    hosts <- vapply(res, function(z) as.character(z$cluster[1L]), "")
    j <- which(hosts == name)
    if (!length(j)) return(NULL)
    res[[j[1L]]]
}


test_that("findNodeMarkers returns markers and a per-node summary", {
    skip_if_not_installed("scran")
    fx <- .nm_gobject()
    tree <- .nm_tree(fx$gobject)
    splits <- getDendrogramSplits(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, show_dend = FALSE,
        verbose = FALSE)

    res <- findNodeMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, splits = splits,
        method = "scran", verbose = FALSE)

    expect_named(res, c("markers", "nodes"))
    expect_true(all(c("nodeID", "side", "feats") %in% names(res$markers)))
    expect_true(all(c("nodeID", "node_h", "left", "right", "n_strong") %in%
        names(res$nodes)))
    expect_type(res$nodes$nodeID, "integer")
    expect_type(res$nodes$node_h, "double")

    # one row per internal node, and every merge row represented exactly once
    expect_identical(nrow(res$nodes), nrow(tree$merge))
    expect_setequal(res$nodes$nodeID, seq_len(nrow(tree$merge)))
    expect_false(anyDuplicated(res$nodes$nodeID) > 0L)

    # side is only ever left or right
    expect_setequal(unique(res$markers$side), c("left", "right"))

    # the two sides of each node are disjoint and partition that subtree
    for (i in seq_len(nrow(splits))) {
        expect_length(intersect(splits$tree_1[[i]], splits$tree_2[[i]]), 0L)
    }
})


test_that("pooled scran node markers equal the delegated findMarkers call", {
    skip_if_not_installed("scran")
    fx <- .nm_gobject()
    tree <- .nm_tree(fx$gobject)
    splits <- getDendrogramSplits(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, show_dend = FALSE,
        verbose = FALSE)
    pooled <- findNodeMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, splits = splits,
        method = "scran", verbose = FALSE)

    for (i in seq_len(nrow(splits))) {
        L <- splits$tree_1[[i]]
        R <- splits$tree_2[[i]]
        nid <- splits$nodeID[i]
        deleg <- findMarkers(fx$gobject, cluster_column = "clus",
            expression_values = "raw", method = "scran",
            group_1 = L, group_2 = R, verbose = FALSE)
        # `findMarkers()` names the group by pasting it unsorted
        host <- .nm_host(deleg, paste0(L, collapse = "_"))
        expect_false(is.null(host))

        lf <- grep("^logFC", names(host), value = TRUE)[1L]
        got <- pooled$markers[
            pooled$markers$nodeID == nid & pooled$markers$side == "left", ]
        cmp <- merge(
            data.table::data.table(feats = got$feats, a = got$logFC),
            data.table::data.table(feats = host$feats, d = host[[lf]]),
            by = "feats")
        expect_gt(nrow(cmp), 0L)
        expect_equal(cmp$a, cmp$d, tolerance = 1e-10)
    }
})


test_that("pooled gini node markers equal the delegated findMarkers call", {
    fx <- .nm_gobject()
    tree <- .nm_tree(fx$gobject)
    splits <- getDendrogramSplits(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, show_dend = FALSE,
        verbose = FALSE)
    pooled <- findNodeMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, splits = splits,
        method = "gini")

    for (i in seq_len(nrow(splits))) {
        L <- splits$tree_1[[i]]
        R <- splits$tree_2[[i]]
        nid <- splits$nodeID[i]
        lname <- paste0(L, collapse = "_")
        deleg <- data.table::as.data.table(findMarkers(fx$gobject,
            cluster_column = "clus", expression_values = "raw",
            method = "gini", group_1 = L, group_2 = R))

        got <- pooled$markers[pooled$markers$nodeID == nid, ]
        # same number of surviving rows: this is what the parameter defaults
        # govern, and what regressed when they were inherited from markersParam
        expect_identical(nrow(got), nrow(deleg))

        # gini keys rows by (feats, group); map the delegated group label to a
        # side so the join is on the same key the pooled table uses
        d <- data.table::data.table(
            feats = deleg$feats,
            side = ifelse(as.character(deleg$cluster) == lname,
                "left", "right"),
            d = deleg$comb_score)
        cmp <- merge(
            data.table::data.table(feats = got$feats, side = got$side,
                a = got$comb_score),
            d, by = c("feats", "side"))
        expect_identical(nrow(cmp), nrow(got))
        # `comb_score` is NaN for a feature absent from both sides; both paths
        # must produce NaN for the same features, not merely similar numbers
        expect_identical(is.na(cmp$a), is.na(cmp$d))
        expect_equal(cmp$a[!is.na(cmp$a)], cmp$d[!is.na(cmp$d)],
            tolerance = 1e-10)
    }
})


test_that("gini node markers use findMarkers defaults, not markersParam's", {
    fx <- .nm_gobject()
    tree <- .nm_tree(fx$gobject)

    # the defaults on the signature
    f <- formals(findNodeMarkers)
    expect_identical(eval(f$min_expression), 0.5)
    expect_identical(eval(f$min_detection), 0.5)
    expect_identical(eval(f$min_feats), 4)

    # and they must actually reach the statistic: a looser threshold has to
    # return at least as many rows as the default one
    tight <- findNodeMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, method = "gini")
    loose <- findNodeMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, method = "gini",
        min_expression = 0.2, min_detection = 0.2, min_feats = 5)
    expect_gte(nrow(loose$markers), nrow(tight$markers))
})


test_that("node markers match scran::findMarkers on a relabelled matrix", {
    skip_if_not_installed("scran")
    fx <- .nm_gobject()
    tree <- .nm_tree(fx$gobject)
    splits <- getDendrogramSplits(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, show_dend = FALSE,
        verbose = FALSE)
    pooled <- findNodeMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, splits = splits,
        method = "scran", verbose = FALSE)

    # ground truth from outside the suite: relabel the cells of one node and
    # hand the matrix straight to scran, sharing no code with the implementation
    i <- which.max(splits$node_h)
    L <- splits$tree_1[[i]]
    R <- splits$tree_2[[i]]
    keep <- fx$clus %in% c(L, R)
    lab <- ifelse(fx$clus[keep] %in% L, "L", "R")
    ref <- scran::findMarkers(as.matrix(fx$mat[, keep, drop = FALSE]),
        groups = lab)
    rl <- data.table::as.data.table(ref[["L"]])
    rl[, "feats" := rownames(ref[["L"]])]

    got <- pooled$markers[
        pooled$markers$nodeID == splits$nodeID[i] &
            pooled$markers$side == "left", ]
    cmp <- merge(
        data.table::data.table(feats = got$feats, a = got$p.value),
        data.table::data.table(feats = rl$feats, d = rl$p.value),
        by = "feats")
    expect_gt(nrow(cmp), 0L)
    expect_equal(cmp$a, cmp$d, tolerance = 1e-8)
})


test_that("findNodeMarkers builds its own tree and splits when not given them", {
    skip_if_not_installed("scran")
    fx <- .nm_gobject()

    # `cor` / `distance` must reach the tree this builds for itself, or the
    # auto path silently uses a different linkage from the explicit one
    auto <- findNodeMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", cor = "pearson", distance = "average",
        method = "scran", verbose = FALSE)
    tree <- .nm_tree(fx$gobject)
    given <- findNodeMarkers(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, method = "scran",
        verbose = FALSE)

    expect_identical(auto$nodes$nodeID, given$nodes$nodeID)
    expect_equal(auto$nodes$node_h, given$nodes$node_h, tolerance = 1e-12)
})


test_that("findNodeMarkers rejects a tree that is not an hclust", {
    fx <- .nm_gobject()
    expect_error(
        findNodeMarkers(fx$gobject, cluster_column = "clus",
            expression_values = "raw", tree = list(a = 1), method = "scran"),
        "must be an `hclust`"
    )
})
