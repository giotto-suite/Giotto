# Building the annotation query for a cluster tree, and reading the answer back
# onto the object at several granularities.
#
# The load-bearing case is a *partially labelled* tree. A real annotator does
# not name every internal node, and the resolution order in
# `annotateClusterTree()` exists so a column still comes back with no holes in
# it. Two bugs here reached a rendered document before being caught, both the
# same shape: `[[` on a named vector errors for an absent name rather than
# returning NULL, so the `%null%` fallback never ran.

.ct_gobject <- function(n_genes = 40L, n_cells = 90L, n_clus = 6L, seed = 5L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.7,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    clus <- as.character(rep(seq_len(n_clus), length.out = n_cells))
    for (k in seq_len(n_clus)) {
        gi <- ((k - 1L) * 4L + 1L):(k * 4L)
        m[gi, clus == as.character(k)] <- m[gi, clus == as.character(k)] + 20
    }
    g <- GiottoClass::createGiottoObject(expression = m)
    g <- GiottoClass::addCellMetadata(g,
        new_metadata = data.frame(cell_ID = colnames(m), clus = clus))
    list(gobject = g, mat = m, clus = clus, n_clus = n_clus)
}

.ct_tree <- function(g) {
    calculateClusterTree(g, cluster_column = "clus",
        expression_values = "raw", cor = "pearson", distance = "average")
}

.ct_labels <- function(tree, nodes = TRUE) {
    out <- list(clusters = stats::setNames(
        paste0("type_", tree$labels), tree$labels
    ))
    if (isTRUE(nodes)) {
        # deliberately partial: only the deepest few nodes are named, which is
        # what a real answer looks like
        ids <- as.character(utils::tail(seq_len(nrow(tree$merge)), 2L))
        out$nodes <- stats::setNames(paste0("clade_", ids), ids)
    }
    out
}


test_that("writeClusterTreeQuery renders the context and all four sections", {
    skip_if_not_installed("scran")
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)

    q <- writeClusterTreeQuery(fx$gobject,
        cluster_column = "clus", expression_values = "raw", tree = tree,
        context = list(tissue = "test tissue", disease = "none",
            assay = "made up")
    )

    expect_type(q, "character")
    expect_true(any(grepl("^tissue: test tissue$", q)))
    expect_true(any(grepl("^disease: none$", q)))
    expect_true(any(grepl("^assay: made up$", q)))

    for (sec in c("## 1. Cluster tree", "## 2. Markers at each split",
                  "## 3. Per-cluster markers", "## 4. Required output")) {
        expect_true(any(startsWith(q, sec)), info = sec)
    }

    # every cluster and every node has to appear, or the answer cannot cover
    # the tree it was asked about
    for (cl in tree$labels) {
        expect_true(any(grepl(paste0("^### Cluster ", cl, " "), q)), info = cl)
    }
    expect_equal(sum(startsWith(q, "### Node ")), nrow(tree$merge))

    # a node label must name the clade, not describe the split -- this is what
    # makes the answer cuttable afterwards
    expect_true(any(grepl("names the group of clusters BENEATH", q)))
})


test_that("writeClusterTreeQuery returns the text it writes", {
    skip_if_not_installed("scran")
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)
    f <- tempfile(fileext = ".txt")

    q <- writeClusterTreeQuery(fx$gobject, cluster_column = "clus",
        expression_values = "raw", tree = tree, file_name = f)

    expect_true(file.exists(f))
    expect_identical(readLines(f), q)
})


test_that("annotateClusterTree writes one column per requested level", {
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)
    labs <- .ct_labels(tree)

    g <- annotateClusterTree(fx$gobject, tree = tree, labels = labs,
        cluster_column = "clus", k = c(2L, 3L, fx$n_clus))
    cd <- GiottoClass::pDataDT(g)

    for (cc in c("cell_types_k2", "cell_types_k3",
                 paste0("cell_types_k", fx$n_clus))) {
        expect_true(cc %in% colnames(cd), info = cc)
        # the fallback chain exists so a partially labelled tree still fills
        # every cell; a hole here means it did not
        expect_false(anyNA(cd[[cc]]), info = cc)
    }
    expect_length(unique(cd$cell_types_k2), 2L)
    expect_length(unique(cd$cell_types_k3), 3L)
})


test_that("the finest cut reproduces the leaf labels exactly", {
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)
    labs <- .ct_labels(tree)

    g <- annotateClusterTree(fx$gobject, tree = tree, labels = labs,
        cluster_column = "clus", k = fx$n_clus)
    cd <- GiottoClass::pDataDT(g)

    # no internal node spans a single leaf, so at the finest cut every group is
    # a singleton. If the ancestor search ran first these would come back as
    # clade names and the finest level would be coarser than its own input.
    expect_identical(
        as.character(cd[[paste0("cell_types_k", fx$n_clus)]]),
        unname(labs$clusters[as.character(cd$clus)])
    )
})


test_that("annotateClusterTree fills every level with no node labels at all", {
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)
    labs <- .ct_labels(tree, nodes = FALSE)

    g <- annotateClusterTree(fx$gobject, tree = tree, labels = labs,
        cluster_column = "clus", k = c(2L, 4L))
    cd <- GiottoClass::pDataDT(g)

    expect_false(anyNA(cd$cell_types_k2))
    expect_false(anyNA(cd$cell_types_k4))
})


test_that("a node id absent from the labels does not error", {
    # the bug this pins: `[[` on a named vector throws "subscript out of
    # bounds" for a missing name, so guarding with `%null%` afterwards is too
    # late. Every k must survive a sparsely labelled node set.
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)
    labs <- list(
        clusters = stats::setNames(paste0("type_", tree$labels), tree$labels),
        nodes = stats::setNames("only_one", "1")
    )
    for (kk in seq_len(fx$n_clus)) {
        expect_error(
            annotateClusterTree(fx$gobject, tree = tree, labels = labs,
                cluster_column = "clus", k = kk),
            NA
        )
    }
})


test_that("annotateClusterTree takes h as well as k", {
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)
    labs <- .ct_labels(tree)
    hh <- stats::median(tree$height)

    g <- annotateClusterTree(fx$gobject, tree = tree, labels = labs,
        cluster_column = "clus", h = hh)
    cd <- GiottoClass::pDataDT(g)

    cc <- grep("^cell_types_h", colnames(cd), value = TRUE)
    expect_length(cc, 1L)
    expect_false(anyNA(cd[[cc]]))
})


test_that("labels are accepted as a table as well as a named vector", {
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)

    as_vec <- list(clusters = stats::setNames(
        paste0("type_", tree$labels), tree$labels))
    as_tab <- list(clusters = data.frame(
        cluster = tree$labels, cell_type = paste0("type_", tree$labels)))

    a <- GiottoClass::pDataDT(annotateClusterTree(fx$gobject, tree = tree,
        labels = as_vec, cluster_column = "clus", k = 3L))
    b <- GiottoClass::pDataDT(annotateClusterTree(fx$gobject, tree = tree,
        labels = as_tab, cluster_column = "clus", k = 3L))
    expect_identical(a$cell_types_k3, b$cell_types_k3)
})


test_that("granularity arguments are refused in Giotto's own terms", {
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)
    labs <- .ct_labels(tree)

    expect_error(
        annotateClusterTree(fx$gobject, tree = tree, labels = labs,
            cluster_column = "clus", k = fx$n_clus + 5L),
        "must be between 1 and the number of clusters"
    )
    expect_error(
        annotateClusterTree(fx$gobject, tree = tree, labels = labs,
            cluster_column = "clus", k = NA_integer_),
        "must be numeric and not NA"
    )
    expect_error(
        annotateClusterTree(fx$gobject, tree = tree, labels = labs,
            cluster_column = "clus", h = -1),
        "must not be negative"
    )
    expect_error(
        annotateClusterTree(fx$gobject, tree = list(1), labels = labs,
            cluster_column = "clus", k = 2L),
        "must be an `hclust`"
    )
})


test_that("missing cluster labels are named in the error", {
    fx <- .ct_gobject()
    tree <- .ct_tree(fx$gobject)
    labs <- list(clusters = stats::setNames("type_1", "1"))

    expect_error(
        annotateClusterTree(fx$gobject, tree = tree, labels = labs,
            cluster_column = "clus", k = 2L),
        "is missing"
    )
})
