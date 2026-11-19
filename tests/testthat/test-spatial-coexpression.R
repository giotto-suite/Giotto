.sc_g <- function() {
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    feats <- spatialAutoCor(g, method = "binspect",
        spatial_network_name = "spatial_network",
        return_gobject = FALSE)$feats[1:60]
    list(g = g, feats = feats)
}

test_that("the result is the four building blocks composed", {
    o <- .sc_g()
    res <- spatialCoexpression(o$g, feats = o$feats,
        spatial_network_name = "spatial_network", smoothing = 0.6,
        n_modules = 4, return_gobject = FALSE)

    x <- getExpression(o$g, values = "normalized", output = "matrix")[o$feats, ]
    net <- getSpatialNetwork(o$g, name = "spatial_network",
        output = "networkDT")
    sm <- processData(x, spatialSmoothParam(smoothing = 0.6), network = net)
    cm <- analyzeData(sm, featCorParam())
    expect_equal(res$cor, cm)
    hc <- stats::hclust(stats::as.dist(1 - cm), method = "ward.D")
    expect_identical(res$modules, stats::cutree(hc, k = 4))
    expect_identical(sort(unique(res$summary$module)), 1:4)
    expect_equal(sum(res$summary$size), length(o$feats))
})

test_that("smoothing = 0 is plain co-expression", {
    o <- .sc_g()
    res <- spatialCoexpression(o$g, feats = o$feats[1:20],
        spatial_network_name = "spatial_network", smoothing = 0,
        n_modules = 3, return_gobject = FALSE)
    x <- getExpression(o$g, values = "normalized", output = "matrix")
    expect_equal(res$cor, cor(t(as.matrix(x[o$feats[1:20], ]))))
})

test_that("the giotto object gets modules per feature and scores per cell", {
    o <- .sc_g()
    g <- spatialCoexpression(o$g, feats = o$feats,
        spatial_network_name = "spatial_network", n_modules = 4)
    fm <- fDataDT(g)
    expect_true(all(c("spatial_module", "spatial_module_cor") %in% names(fm)))
    expect_identical(sort(fm[!is.na(spatial_module), feat_ID]), sort(o$feats))
    enr <- getSpatialEnrichment(g, name = "spatial_module",
        output = "data.table")
    expect_identical(setdiff(names(enr), "cell_ID"), as.character(1:4))
    expect_identical(nrow(enr), length(spatIDs(g)))
})

test_that("heatmapCoexpression draws the modules", {
    o <- .sc_g()
    res <- spatialCoexpression(o$g, feats = o$feats,
        spatial_network_name = "spatial_network", n_modules = 4,
        return_gobject = FALSE)
    h <- heatmapCoexpression(o$g, res, show_plot = FALSE, return_plot = TRUE,
        save_plot = FALSE)
    expect_s4_class(h, "Heatmap")
    expect_error(heatmapCoexpression(o$g, list(a = 1)), "spatialCoexpression")
})
