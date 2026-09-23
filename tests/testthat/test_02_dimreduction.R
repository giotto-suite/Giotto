

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

# umap reproducibility ####

# There was no test anywhere in the suite asserting that runUMAP() reproduces
# across two calls, which is why it did not. It carried set_seed = TRUE and
# still differed run to run, because uwot ran its own approximate search whose
# HNSW index build races on insertion order -- a C++ thread interleaving the R
# RNG cannot reach. `nn_engine = "giotto"` supplies the graph instead.

test_that("runUMAP is reproducible at default threads", {
    skip_if_not_installed("uwot")
    skip_if_not_installed("RcppHNSW")
    rlang::local_options(lifecycle_verbosity = "quiet")

    a <- suppressWarnings(runUMAP(g, name = "umap_rep_a", n_neighbors = 15))
    b <- suppressWarnings(runUMAP(g, name = "umap_rep_b", n_neighbors = 15))

    ca <- getDimReduction(a, reduction_method = "umap",
                          name = "umap_rep_a", output = "dimObj")[]
    cb <- getDimReduction(b, reduction_method = "umap",
                          name = "umap_rep_b", output = "dimObj")[]
    expect_identical(ca, cb)
})

test_that("nn_engine defaults to the graph Giotto builds itself", {
    expect_identical(eval(formals(runUMAP)$nn_engine)[1], "giotto")
    # every uwot backend is reachable through the same argument, rather than
    # annoy being privileged because it is the one that got measured
    expect_setequal(eval(formals(runUMAP)$nn_engine),
        c("giotto", "uwot", "fnn", "annoy", "hnsw", "nndescent"))
})

test_that("a stored kNN is chosen by name, not by being listed first", {
    skip_if_not_installed("uwot")
    skip_if_not_installed("RcppHNSW")
    rlang::local_options(lifecycle_verbosity = "quiet")

    # Two kNN networks whose NEIGHBOURS genuinely differ. Varying only k
    # would not do it: the nearest 9 of a 12-neighbour graph are the same
    # nine as the nearest 9 of a 25-neighbour graph, so both would embed
    # identically and the test would pass for the wrong reason.
    gg <- createNearestNetwork(g, type = "kNN", dim_reduction_to_use = "pca",
        dimensions_to_use = 1:10, k = 12, name = "kNN.pca", verbose = FALSE)
    gg <- createNearestNetwork(gg, type = "kNN", dim_reduction_to_use = "pca",
        dimensions_to_use = 1:2, k = 12, name = "kNN.flat", verbose = FALSE)
    expect_true(all(c("kNN.pca", "kNN.flat") %in%
                    list_nearest_networks(gg)$name))

    # the named network is the one reported as used
    expect_message(
        runUMAP(gg, name = "u_msg", n_neighbors = 10,
                network_name = "kNN.flat", verbose = TRUE),
        "kNN.flat"
    )

    deep <- suppressWarnings(runUMAP(gg, name = "u_deep",
        n_neighbors = 10, network_name = "kNN.pca"))
    flat <- suppressWarnings(runUMAP(gg, name = "u_flat",
        n_neighbors = 10, network_name = "kNN.flat"))
    cd <- getDimReduction(deep, reduction_method = "umap",
                          name = "u_deep", output = "dimObj")[]
    cf <- getDimReduction(flat, reduction_method = "umap",
                          name = "u_flat", output = "dimObj")[]
    # different neighbourhoods, so different embeddings -- were the name
    # ignored these would agree and the argument would be decorative
    expect_false(identical(cd, cf))

    # and each is reproducible on its own named graph
    again <- suppressWarnings(runUMAP(gg, name = "u_deep2",
        n_neighbors = 10, network_name = "kNN.pca"))
    expect_identical(cd, getDimReduction(again, reduction_method = "umap",
        name = "u_deep2", output = "dimObj")[])
})

test_that("an explicitly named network that is absent is an error", {
    skip_if_not_installed("uwot")
    rlang::local_options(lifecycle_verbosity = "quiet")
    # falling back quietly would embed on a graph the caller did not ask for
    expect_error(
        runUMAP(g, name = "u_missing", n_neighbors = 10,
                network_name = "kNN.does_not_exist"),
        "no kNN network named"
    )
})

test_that("uwot's own backends are reachable through nn_engine", {
    skip_if_not_installed("uwot")
    rlang::local_options(lifecycle_verbosity = "quiet")
    for (be in c("fnn", "annoy")) {
        u <- suppressWarnings(runUMAP(g, name = paste0("u_", be),
            n_neighbors = 10, nn_engine = be))
        cu <- getDimReduction(u, reduction_method = "umap",
                              name = paste0("u_", be), output = "dimObj")[]
        checkmate::expect_matrix(cu, ncols = 2L)
        expect_false(anyNA(cu))
    }
})

test_that("nn_engine = 'annoy' is reproducible too", {
    skip_if_not_installed("uwot")
    rlang::local_options(lifecycle_verbosity = "quiet")

    a <- suppressWarnings(runUMAP(g, name = "umap_ann_a", n_neighbors = 15,
                                  nn_engine = "annoy"))
    b <- suppressWarnings(runUMAP(g, name = "umap_ann_b", n_neighbors = 15,
                                  nn_engine = "annoy"))
    expect_identical(
        getDimReduction(a, reduction_method = "umap",
                        name = "umap_ann_a", output = "dimObj")[],
        getDimReduction(b, reduction_method = "umap",
                        name = "umap_ann_b", output = "dimObj")[]
    )
})

test_that("nn_engine = 'giotto' supplies the graph and is reproducible", {
    skip_if_not_installed("uwot")
    skip_if_not_installed("RcppHNSW")
    rlang::local_options(lifecycle_verbosity = "quiet")

    a <- suppressWarnings(runUMAP(g, name = "umap_gio_a", n_neighbors = 15,
                                  nn_engine = "giotto"))
    b <- suppressWarnings(runUMAP(g, name = "umap_gio_b", n_neighbors = 15,
                                  nn_engine = "giotto"))
    ca <- getDimReduction(a, reduction_method = "umap",
                          name = "umap_gio_a", output = "dimObj")[]
    cb <- getDimReduction(b, reduction_method = "umap",
                          name = "umap_gio_b", output = "dimObj")[]
    expect_identical(ca, cb)
})

test_that("runUMAP refuses a caller graph while building its own", {
    skip_if_not_installed("uwot")
    rlang::local_options(lifecycle_verbosity = "quiet")
    # silently dropping the caller's graph would be the worse failure: they
    # would get an embedding of a different graph than the one they passed
    expect_error(
        runUMAP(g, nn_method = list(idx = 1, dist = 1)),
        "`nn_method` was supplied"
    )
    expect_error(
        runUMAP(g, nn_method = "hnsw", nn_engine = "giotto"),
        "`nn_method` was supplied"
    )
})

test_that("nn_engine = 'uwot' still reaches uwot's own search", {
    skip_if_not_installed("uwot")
    rlang::local_options(lifecycle_verbosity = "quiet")
    u <- suppressWarnings(
        runUMAP(g, name = "umap_uwot", n_neighbors = 15, nn_engine = "uwot")
    )
    cu <- getDimReduction(u, reduction_method = "umap",
                          name = "umap_uwot", output = "dimObj")[]
    checkmate::expect_matrix(cu, ncols = 2L)
    expect_false(anyNA(cu))
})
