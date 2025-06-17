

# DATA TO USE
# visium mini expression
g <- test_data$vis

# pca ####
test_that("pca is calculated", {
    rlang::local_options(lifecycle_verbosity = "quiet")
    # remove dim reductions
    g@dimension_reduction <- NULL
    g@nn_network <- NULL

    g <- suppressWarnings(runPCA(g))
    x <- getDimReduction(g)

    checkmate::expect_class(x, "dimObj")
    checkmate::expect_numeric(x$eigenvalues)
    checkmate::expect_matrix(x$loadings)
})

test_that("projection pca is calculated", {
    rlang::local_options(lifecycle_verbosity = "quiet")
    # remove dim reductions
    g@dimension_reduction <- NULL
    g@nn_network <- NULL

    g <- suppressWarnings(runPCAprojection(g))
    x <- getDimReduction(g, name = "pca.projection")

    checkmate::expect_class(x, "dimObj")
    checkmate::expect_numeric(x$eigenvalues)
    checkmate::expect_matrix(x$loadings)
})

# UMAP ####
test_that("umap is calculated", {
    rlang::local_options(lifecycle_verbosity = "quiet")
    g <- setDimReduction(g, NULL,
        spat_unit = "cell",
        feat_type = "rna",
        reduction = "cells",
        reduction_method = "umap",
        name = "umap"
    )
    expect_false("umap" %in% list_dim_reductions(g)$name)

    g <- runUMAP(g)

    expect_true("umap" %in% list_dim_reductions(g)$name)

    u <- getDimReduction(g,
        spat_unit = "cell",
        feat_type = "rna",
        reduction = "cells",
        reduction_method = "umap",
        name = "umap"
    )

    checkmate::expect_class(u, "dimObj")
})
