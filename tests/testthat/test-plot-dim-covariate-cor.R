.pdcc_g <- function() {
    GiottoData::loadGiottoMini("visium", verbose = FALSE)
}

test_that("one covariate gives the absolute correlation of each component", {
    g <- .pdcc_g()
    p <- plotDimCovariateCor(g, covariates = "total_expr",
        dimensions_to_use = 10, show_plot = FALSE, return_plot = TRUE,
        save_plot = FALSE)
    pcs <- getDimReduction(g, reduction_method = "pca", output = "matrix")
    cm <- pDataDT(g)
    x <- log10(pmax(cm$total_expr[match(rownames(pcs), cm$cell_ID)], 1))
    expect_equal(p$data$abs_cor, abs(as.numeric(stats::cor(pcs[, 1:10], x))))
    expect_match(p$labels$subtitle, "^PC1 \\|r\\| = ")
})

test_that("transform is given once or per covariate", {
    g <- .pdcc_g()
    p <- plotDimCovariateCor(g, covariates = c("total_expr", "nr_feats"),
        transform = c("log10", "none"), dimensions_to_use = 5,
        show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
    pcs <- getDimReduction(g, reduction_method = "pca", output = "matrix")
    cm <- pDataDT(g)
    nf <- cm$nr_feats[match(rownames(pcs), cm$cell_ID)]
    expect_equal(p$data[p$data$covariate == "nr_feats", ]$abs_cor,
        abs(as.numeric(stats::cor(pcs[, 1:5], nf))))
    expect_error(plotDimCovariateCor(g, covariates = c("total_expr",
        "nr_feats"), transform = c("log10", "none", "none")),
        "one per covariate")
})

test_that("covariates can be a matrix keyed by cell ID", {
    g <- .pdcc_g()
    sl <- getSpatialLocations(g, output = "data.table")
    xy <- as.matrix(sl[, c("sdimx", "sdimy")])
    rownames(xy) <- sl$cell_ID
    p <- plotDimCovariateCor(g, covariates = xy, transform = "none",
        dimensions_to_use = 5, show_plot = FALSE, return_plot = TRUE,
        save_plot = FALSE)
    expect_setequal(unique(as.character(p$data$covariate)),
        c("sdimx", "sdimy"))
})

test_that("unusable inputs stop with a message", {
    g <- .pdcc_g()
    expect_error(plotDimCovariateCor(g, dim_reduction_to_use = "umap"),
        "not a linear reduction")
    expect_error(plotDimCovariateCor(g, covariates = "not_a_column"),
        "not in cell metadata")
    expect_error(plotDimCovariateCor(g, covariates = "cell_ID"),
        "not numeric")
})
