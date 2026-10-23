.sck_g <- function() {
    GiottoData::loadGiottoMini("visium", verbose = FALSE)
}

test_that("wss on the PCA is the k-means objective at each k", {
    g <- .sck_g()
    tb <- sweepClusterK(g, k_range = 2:6, score = "wss", seed = 1,
        output = "data.table")
    pcs <- getDimReduction(g, reduction_method = "pca", output = "matrix")[,
        1:10]
    set.seed(1)
    wss <- vapply(2:6, function(k) stats::kmeans(pcs, centers = k, nstart = 5,
        iter.max = 50)$tot.withinss, numeric(1))
    expect_equal(tb$value, wss)
    expect_identical(names(tb), c("k", "value", "se", "score", "best"))
    expect_identical(sum(tb$best), 1L)
})

test_that("a co-expression result is swept along its own tree", {
    g <- .sck_g()
    coexp <- spatialCoexpression(g, feats = head(rownames(g), 60),
        spatial_network_name = "spatial_network", return_gobject = FALSE)
    tb <- sweepClusterK(coexp, k_range = 2:8, output = "data.table")
    d <- stats::as.dist(1 - coexp$cor)
    sil <- vapply(2:8, function(k) mean(cluster::silhouette(
        stats::cutree(coexp$hclust, k), d)[, 3L]), numeric(1))
    expect_equal(tb$value, sil)
    expect_identical(unique(tb$score), "silhouette")

    p <- sweepClusterK(coexp, k_range = 2:8, show_plot = FALSE,
        return_plot = TRUE)
    expect_identical(p$labels$x, "number of modules (n_modules)")
})

test_that("a spatial enrichment is swept for niches", {
    g <- .sck_g()
    g <- calculateLabelProportions(g, labels = "leiden_clus",
        group_method = "spatialnetwork",
        spatial_network_name = "spatial_network", name = "niche_props",
        output = "gobject")
    tb <- sweepClusterK(g, values_to_use = "spatial_enrichment",
        spat_enr_name = "niche_props", k_range = 2:6, score = "wss", seed = 1,
        output = "data.table")
    expect_identical(tb$k, 2:6)
})

test_that("inputs that cannot be scored stop with a message", {
    g <- .sck_g()
    coexp <- spatialCoexpression(g, feats = head(rownames(g), 60),
        spatial_network_name = "spatial_network", return_gobject = FALSE)
    expect_error(sweepClusterK(coexp, k_range = 2:8, score = "wss"),
        "only 'silhouette' applies")
    expect_error(sweepClusterK(g, dim_reduction_to_use = "umap"),
        "not a linear reduction")
    expect_error(sweepClusterK(g, k_range = 2:3), "k_range")
})

test_that("a setting the input does not use warns only when it was set", {
    g <- .sck_g()
    coexp <- spatialCoexpression(g, feats = head(rownames(g), 60),
        spatial_network_name = "spatial_network", return_gobject = FALSE)
    expect_no_warning(sweepClusterK(coexp, k_range = 2:8, show_plot = FALSE))
    expect_warning(sweepClusterK(coexp, k_range = 2:8,
        hclust_method = "average", show_plot = FALSE), "own tree")
    expect_warning(sweepClusterK(g, k_range = 2:6, score = "wss", seed = 1,
        spat_enr_name = "niche_props", show_plot = FALSE),
        "values_to_use = \"spatial_enrichment\"")
})
