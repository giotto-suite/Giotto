
test_that("leiden (python) clustering works", {
    g <- test_data$vis
    test <- doLeidenCluster(g,
        resolution = 0.1,
        name = "test_col"
    )
    cmeta <- pDataDT(test)
    expect_true("test_col" %in% colnames(cmeta))
    expect_true(cmeta$test_col[200] == 1)
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
