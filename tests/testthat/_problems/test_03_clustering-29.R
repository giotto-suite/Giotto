# Extracted from test_03_clustering.R:29

# test -------------------------------------------------------------------------
g <- test_data$vis
test <- doLouvainCluster(g,version = "community",
        resolution = 0.1,
        name = "test_col"
    )
