test_that("leiden through clusterCells and subclustering works", {
    g <- test_data$vis
    test <- clusterCells(g, cluster_method = "leiden",
        pyth_leid_resolution = 0.1, name = "test_col")
    expect_true("test_col" %in% colnames(pDataDT(test)))
    dt <- doLeidenCluster(g, resolution = 0.1, n_iterations = 1000,
        name = "test_col",
        return_gobject = FALSE)
    expect_named(dt, c("cell_ID", "test_col"))
    cm <- pDataDT(test)
    expect_identical(dt$test_col[match(cm$cell_ID, dt$cell_ID)], cm$test_col)
    sub <- doLeidenSubCluster(g, cluster_column = "leiden_clus",
        selected_clusters = c(1, 2), name = "sub_col")
    sub_ids <- unique(pDataDT(sub)$sub_col)
    expect_true(any(startsWith(sub_ids, "1.")))
    expect_warning(
        clusterCells(g, cluster_method = "leiden", python_path = "x"),
        class = "lifecycle_warning_deprecated"
    )
})

test_that("leiden (igraph) clustering works", {
    g <- test_data$vis
    test <- doLeidenClusterIgraph(g,
        resolution = 0.1,
        name = "test_col"
    )
    cmeta <- pDataDT(test)
    expect_true("test_col" %in% colnames(cmeta))
    expect_true(cmeta$test_col[200] == 1)
})

test_that("louvain (igraph) clustering works", {
    g <- test_data$vis
    test <- doLouvainCluster(g,
        resolution = 0.5,
        name = "test_col"
    )
    cmeta <- pDataDT(test)
    expect_true("test_col" %in% colnames(cmeta))
    expect_gt(length(unique(cmeta$test_col)), 1L)
    # the default seed makes the clusters reproducible
    again <- pDataDT(doLouvainCluster(g, resolution = 0.5, name = "test_col"))
    expect_identical(again$test_col, cmeta$test_col)
    # higher resolution, more clusters
    fine <- pDataDT(doLouvainCluster(g, resolution = 2, name = "test_col"))
    expect_gt(length(unique(fine$test_col)), length(unique(cmeta$test_col)))
    expect_error(clusterParam("louvain_community"))
})

# multinet testing skipped
