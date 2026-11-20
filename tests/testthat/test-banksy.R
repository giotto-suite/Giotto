skip_if_not_installed("Banksy")
skip_if_not_installed("SpatialExperiment")

.banksy_fixture <- function(backend = NULL) {
    withr::local_seed(3)
    n <- 400
    cells <- sprintf("c_%03d", seq_len(n))
    x <- stats::runif(n, 0, 100)
    y <- stats::runif(n, 0, 100)
    region <- as.integer(x > 50) + 1L
    emb <- matrix(stats::rnorm(n * 6), n, 6,
        dimnames = list(cells, paste0("Dim.", 1:6)))
    emb[, 1] <- emb[, 1] + 4 * region
    m <- Matrix::Matrix(matrix(stats::rpois(5 * n, 3), 5, n,
        dimnames = list(paste0("g", 1:5), cells)), sparse = TRUE)
    g <- createGiottoObject(expression = m,
        spatial_locs = data.frame(sdimx = x, sdimy = y, cell_ID = cells),
        backend = backend, verbose = FALSE)
    g <- setDimReduction(g, createDimObj(coordinates = emb, name = "pca",
        spat_unit = "cell", feat_type = "rna", method = "pca"),
        spat_unit = "cell", feat_type = "rna", reduction = "cells",
        reduction_method = "pca", name = "pca", verbose = FALSE)
    list(g = g, emb = emb, x = x, y = y, region = region)
}

test_that("runBanksy stores Banksy's PCs and, with cluster_method = 'banksy', its domains", {
    fx <- .banksy_fixture()
    g <- runBanksy(fx$g, npcs = 5, resolution = 0.3, cluster_method = "banksy",
        verbose = FALSE)
    md <- pDataDT(g)
    expect_true("banksy" %in% colnames(md))
    expect_false(anyNA(md$banksy))
    pcs <- getDimReduction(g, reduction_method = "banksy", name = "banksy",
        output = "matrix")
    expect_identical(dim(pcs), c(400L, 5L))
    expect_identical(rownames(pcs), md$cell_ID)
    tab <- table(md$banksy, fx$region[match(md$cell_ID, rownames(fx$emb))])
    expect_gt(sum(apply(tab, 1, max)) / sum(tab), 0.9)

    se <- SpatialExperiment::SpatialExperiment(
        assays = list(embedding = t(fx$emb)),
        spatialCoords = cbind(sdimx = fx$x, sdimy = fx$y))
    se <- Banksy::computeBanksy(se, assay_name = "embedding",
        compute_agf = FALSE, k_geom = 30, verbose = FALSE)
    se <- Banksy::runBanksyPCA(se, use_agf = FALSE, lambda = 0.8, npcs = 5,
        seed = 1234)
    se <- Banksy::clusterBanksy(se, use_agf = FALSE, lambda = 0.8, npcs = 5,
        resolution = 0.3, seed = 1234)
    ref <- SummarizedExperiment::colData(se)
    ref <- as.character(ref[[grep("^clust_", colnames(ref))]])
    expect_identical(md$banksy[match(rownames(fx$emb), md$cell_ID)], ref)

    expect_equal(unname(pcs[rownames(fx$emb), ]), unname(
        SingleCellExperiment::reducedDim(se, "PCA_M0_lam0.8")), tolerance = 1e-12,
        ignore_attr = TRUE)

    agf <- runBanksy(fx$g, k_geom = c(15, 30), compute_agf = TRUE, npcs = 5,
        cluster_method = "banksy", name = "agf", verbose = FALSE)
    expect_identical(dim(getDimReduction(agf, reduction_method = "banksy",
        name = "agf", output = "matrix")), c(400L, 5L))
})

test_that("the default clusters the BANKSY PCs with Giotto's sNN and Leiden", {
    fx <- .banksy_fixture()
    g <- runBanksy(fx$g, npcs = 5, resolution = 0.3, k_neighbors = 30,
        verbose = FALSE)
    md <- pDataDT(g)
    tab <- table(md$banksy, fx$region[match(md$cell_ID, rownames(fx$emb))])
    expect_gt(sum(apply(tab, 1, max)) / sum(tab), 0.9)
    expect_true("banksy" %in% list_nearest_networks(g)$name)
    ref <- createNearestNetwork(g, type = "sNN", dim_reduction_to_use = "banksy",
        dim_reduction_name = "banksy", dimensions_to_use = 1:5, k = 30,
        keep_knn = FALSE, name = "ref", verbose = FALSE)
    ref <- doLeidenCluster(ref, name = "ref", network_name = "ref",
        resolution = 0.3, n_iterations = -1, seed_number = 1234)
    expect_identical(pDataDT(ref)$ref, md$banksy)
})

test_that("runBanksy takes dimensions, a name, and checks its input", {
    fx <- .banksy_fixture()
    g <- runBanksy(fx$g, dimensions_to_use = 1:3, name = "dom", npcs = 4,
        verbose = FALSE)
    expect_true("dom" %in% colnames(pDataDT(g)))
    d <- getDimReduction(g, reduction_method = "banksy", name = "dom",
        output = "dimObj")
    expect_identical(d@misc$input$dimensions_to_use, paste0("Dim.", 1:3))
    expect_error(runBanksy(fx$g, dimensions_to_use = 50:60),
        "\\[runBanksy\\] no stored dimension")
    expect_error(runBanksy(fx$g, lambda = 2))
})

test_that("runBanksy works on a GiottoDisk-backed object", {
    skip_if_not_installed("GiottoDisk")
    vault <- withr::local_tempdir()
    mem <- runBanksy(.banksy_fixture()$g, npcs = 5, verbose = FALSE)
    disk <- runBanksy(.banksy_fixture(file.path(vault, "proj"))$g, npcs = 5,
        verbose = FALSE)
    expect_true(inherits(getExpression(disk)[], "dataStore"))
    expect_identical(pDataDT(disk)$banksy, pDataDT(mem)$banksy)
})
