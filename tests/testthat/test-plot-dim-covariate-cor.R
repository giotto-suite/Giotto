.pdcc_g <- function() {
    GiottoData::loadGiottoMini("visium", verbose = FALSE)
}

test_that("one covariate gives the absolute correlation of each component", {
    g <- .pdcc_g()
    p <- plotDimCovariateCor(g, covariates = "total_expr",
        dimensions_to_use = 1:10, show_plot = FALSE, return_plot = TRUE,
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
        transform = c("log10", "none"), dimensions_to_use = 1:5,
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
        dimensions_to_use = 1:5, show_plot = FALSE, return_plot = TRUE,
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

test_that("dimensions_to_use takes indices and drops those not stored", {
    g <- .pdcc_g()
    pcs <- getDimReduction(g, reduction_method = "pca", output = "matrix")
    p <- plotDimCovariateCor(g, covariates = "total_expr",
        dimensions_to_use = c(2:4, ncol(pcs) + 1L), show_plot = FALSE,
        return_plot = TRUE, save_plot = FALSE)
    expect_identical(p$data$dim, 2:4)
    expect_match(p$labels$subtitle, "^PC2 \\|r\\| = ")
    expect_error(plotDimCovariateCor(g, dimensions_to_use = ncol(pcs) + 1L),
        "none of dimensions_to_use")
})

test_that("a covariate matrix is not log-transformed by default", {
    g <- .pdcc_g()
    sl <- getSpatialLocations(g, output = "data.table")
    xy <- as.matrix(sl[, c("sdimx", "sdimy")])
    rownames(xy) <- sl$cell_ID
    p <- plotDimCovariateCor(g, covariates = xy, dimensions_to_use = 1:5,
        show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
    pcs <- getDimReduction(g, reduction_method = "pca", output = "matrix")
    expect_equal(p$data[p$data$covariate == "sdimx", ]$abs_cor,
        abs(as.numeric(stats::cor(pcs[, 1:5],
            xy[rownames(pcs), "sdimx"]))))
})

test_that("a view correlates only the cells it keeps", {
    g <- .pdcc_g()
    md <- pDataDT(g)
    thr <- stats::median(md$total_expr)
    gv <- subset(g, total_expr > thr, view = "v")
    kept <- md[total_expr > thr]$cell_ID
    args <- list(covariates = "total_expr", dimensions_to_use = 1:5,
        show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
    got <- do.call(plotDimCovariateCor, c(list(gv, view = "v"), args))
    ref <- do.call(plotDimCovariateCor,
        c(list(subsetGiotto(g, cell_ids = kept)), args))
    expect_equal(got$data$abs_cor, ref$data$abs_cor)
    full <- do.call(plotDimCovariateCor, c(list(g), args))
    expect_false(isTRUE(all.equal(got$data$abs_cor, full$data$abs_cor)))
})
