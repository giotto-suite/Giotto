# giottoTree: an `hclust` over clusters that records what it was built from.
# Consumers take their defaults from it, an explicit argument wins with a
# warning when it differs, a plain `hclust` needs explicit arguments, and a
# tree whose leaves do not match the data is refused rather than run through.

.gt_gobject <- function(n_genes = 40L, n_cells = 120L, n_clus = 6L, seed = 3L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.7,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    clus <- as.character(rep(seq_len(n_clus), length.out = n_cells))
    for (k in seq_len(n_clus)) {
        gi <- ((k - 1L) * 5L + 1L):(k * 5L)
        m[gi, clus == as.character(k)] <- m[gi, clus == as.character(k)] + 20
    }
    g <- GiottoClass::createGiottoObject(expression = m)
    GiottoClass::addCellMetadata(g, new_metadata = data.frame(
        cell_ID = colnames(m), clus = clus, clus2 = clus,
        keep = clus %in% c("1", "2", "3")
    ))
}

.gt_tree <- function(g, ...) {
    calculateClusterTree(g, cluster_column = "clus",
        expression_values = "raw", distance = "average", ...)
}


test_that("calculateClusterTree returns a giottoTree that is still an hclust", {
    g <- .gt_gobject()
    tree <- .gt_tree(g)
    expect_s3_class(tree, c("giottoTree", "hclust"), exact = TRUE)
    expect_length(stats::cutree(tree, k = 3), length(tree$labels))
    expect_s3_class(stats::as.dendrogram(tree), "dendrogram")
    expect_identical(attr(tree, "params")$cluster_column, "clus")
})

test_that("as.data.table / as.data.frame give one row per split", {
    g <- .gt_gobject()
    tree <- .gt_tree(g)
    dt <- data.table::as.data.table(tree)
    expect_s3_class(dt, "data.table")
    expect_named(dt, c("nodeID", "node_h", "left", "right"))
    expect_identical(nrow(dt), nrow(tree$merge))
    # root first, and the root separates every leaf
    expect_setequal(c(dt$left[[1]], dt$right[[1]]), tree$labels)
    df <- as.data.frame(tree)
    expect_false(inherits(df, "data.table"))
    expect_identical(df$nodeID, dt$nodeID)
})

test_that("getDendrogramSplits is deprecated and still returns its table", {
    g <- .gt_gobject()
    tree <- .gt_tree(g)
    op <- options(lifecycle_verbosity = "warning")
    on.exit(options(op), add = TRUE)
    expect_warning(
        old <- getDendrogramSplits(g, cluster_column = "clus", tree = tree,
            show_dend = FALSE, verbose = FALSE),
        "deprecated"
    )
    new <- data.table::as.data.table(tree)
    expect_identical(old$nodeID, new$nodeID)
    expect_identical(old$tree_1, new$left)
    expect_identical(old$tree_2, new$right)
})

test_that("consumers default cluster_column from the tree", {
    g <- .gt_gobject()
    tree <- .gt_tree(g)
    labs <- list(clusters = stats::setNames(paste0("t", tree$labels),
        tree$labels))
    a <- annotateClusterTree(g, tree, labs, k = 3)
    b <- annotateClusterTree(g, tree, labs, cluster_column = "clus", k = 3)
    expect_identical(pDataDT(a)$cell_types_k3, pDataDT(b)$cell_types_k3)
})

test_that("an explicit value that differs from the tree's warns, and is used", {
    g <- .gt_gobject()
    tree <- .gt_tree(g)
    labs <- list(clusters = stats::setNames(paste0("t", tree$labels),
        tree$labels))
    expect_warning(
        out <- annotateClusterTree(g, tree, labs, cluster_column = "clus2",
            k = 3, name = "via_clus2"),
        "the tree was built with `cluster_column = \"clus\"`"
    )
    expect_true("via_clus2" %in% colnames(pDataDT(out)))
})

test_that("a plain hclust needs cluster_column", {
    g <- .gt_gobject()
    tree <- .gt_tree(g)
    plain <- tree
    class(plain) <- "hclust"
    attr(plain, "params") <- NULL
    labs <- list(clusters = stats::setNames(paste0("t", tree$labels),
        tree$labels))
    expect_error(annotateClusterTree(g, plain, labs, k = 3),
        "`cluster_column` is needed")
    expect_no_error(annotateClusterTree(g, plain, labs,
        cluster_column = "clus", k = 3))
})

test_that("a tree's view is the default, and another view is flagged", {
    g <- .gt_gobject()
    g <- subset(g, keep == TRUE, view = "v")
    g <- subset(g, clus != "none", view = "all")
    tree <- .gt_tree(g, view = "v")
    expect_setequal(tree$labels, c("1", "2", "3"))

    expect_message(
        inh <- findClusterTreeMarkers(g, tree, method = "gini"),
        "using view \"v\" recorded on the tree"
    )
    same <- findClusterTreeMarkers(g, tree, view = "v", method = "gini")
    expect_identical(inh$nodes, same$nodes)

    # a different view is the caller's call, but it is flagged -- and here
    # it brings back clusters the tree has no leaf for, so it is refused
    expect_error(expect_warning(
        findClusterTreeMarkers(g, tree, view = "all", method = "gini"),
        "the tree was built with `view = \"v\"`"
    ), "clusters with no leaf: 4, 5, 6")
})

test_that("a tree from another clustering is refused", {
    g <- .gt_gobject()
    tree <- .gt_tree(g)
    cx <- getCellMetadata(g, output = "cellMetaObj")
    cx[][clus == "6", "clus" := "7"]
    g2 <- setCellMetadata(g, cx, verbose = FALSE, initialize = FALSE)
    expect_error(findClusterTreeMarkers(g2, tree, method = "gini"),
        "clusters with no leaf: 7.*leaves with no cells: 6")
})

test_that("the tree records which features it was built from", {
    g <- .gt_gobject()
    expect_null(attr(.gt_tree(g), "params")$feats)
    # only features present are recorded, in the order the matrix used them
    want <- c("g3", "g1", "nope", "g7")
    tree <- .gt_tree(g, feats = want)
    expect_setequal(attr(tree, "params")$feats, c("g1", "g3", "g7"))
    expect_identical(attr(tree, "params")$n_feats, 3L)
})

test_that("plot(tree) draws the tree, and what = 'heatmap' its own matrix", {
    skip_if_not_installed("ComplexHeatmap")
    g <- .gt_gobject()
    tree <- .gt_tree(g)
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)

    expect_no_warning(expect_identical(plot(tree, hang = -1), tree))

    hm <- plot(tree, what = "heatmap")
    expect_s4_class(hm, "Heatmap")
    # the matrix the branches were built from, so the two cannot disagree
    expect_identical(hm@matrix,
        attr(tree, "cor_matrix")[tree$labels, tree$labels])
    expect_identical(hm@row_dend_param$obj$merge, tree$merge)
})

test_that("a tree without its matrix says where to go instead", {
    g <- .gt_gobject()
    tree <- .gt_tree(g)
    attr(tree, "cor_matrix") <- NULL
    expect_error(plot(tree, what = "heatmap"), "cluster_custom_order")
})
