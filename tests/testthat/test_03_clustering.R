
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

test_that("leiden clustering works", {
    g <- test_data$vis
    test <- doLouvainCluster(g,version = "community",
        resolution = 0.1,
        name = "test_col"
    )
    cmeta <- pDataDT(test)
    expect_true("test_col" %in% colnames(cmeta))
    expect_true(cmeta$test_col[200] == 25)
})

# multinet testing skipped
