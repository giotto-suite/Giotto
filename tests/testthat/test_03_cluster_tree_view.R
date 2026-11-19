# calculateClusterTree(view = ): the tree is built from the cells that survive
# the view. The reference is the same tree built on a subset object holding
# exactly those cells, with the surviving cells worked out from the metadata
# independently of the view machinery.

g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
md <- pDataDT(g)
CLUS <- "leiden_clus"

.same_tree <- function(a, b) {
    expect_identical(a$labels, b$labels)
    expect_identical(a$merge, b$merge)
    expect_equal(a$height, b$height)
    expect_identical(a$order, b$order)
    expect_equal(attr(a, "cor_matrix"), attr(b, "cor_matrix"))
}

test_that("a view that narrows cells matches a tree built on those cells", {
    thr <- stats::median(md$total_expr)
    gv <- subset(g, total_expr > thr, view = "v")
    kept <- md[total_expr > thr]$cell_ID
    # every cluster survives, so only the means move
    expect_setequal(unique(md[cell_ID %in% kept][[CLUS]]), unique(md[[CLUS]]))

    got <- calculateClusterTree(gv, cluster_column = CLUS, view = "v")
    ref <- calculateClusterTree(subsetGiotto(g, cell_ids = kept),
        cluster_column = CLUS)
    .same_tree(got, ref)
    expect_identical(attr(got, "params")$view, "v")

    # and it is not simply the unnarrowed tree
    full <- calculateClusterTree(g, cluster_column = CLUS)
    expect_false(isTRUE(all.equal(
        attr(got, "cor_matrix"), attr(full, "cor_matrix")
    )))
})

test_that("a cluster with no cells left in the view is not a leaf", {
    keep_clus <- sort(unique(md[[CLUS]]))[1:3]
    gv <- subset(g, leiden_clus %in% keep_clus, view = "v")
    got <- calculateClusterTree(gv, cluster_column = CLUS, view = "v")
    expect_setequal(got$labels, as.character(keep_clus))
})

test_that("without a view the tree is unchanged and records view = NULL", {
    a <- calculateClusterTree(g, cluster_column = CLUS)
    expect_null(attr(a, "params")$view)
    expect_true("view" %in% names(attr(a, "params")))
})
