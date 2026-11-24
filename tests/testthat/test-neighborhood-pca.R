options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
    giotto.no_python_warn = TRUE)

set.seed(1)
n <- 60
locs <- matrix(runif(2 * n), ncol = 2,
    dimnames = list(paste0("c", seq_len(n)), c("sdimx", "sdimy")))
x <- matrix(rpois(10 * n, 3), nrow = 10,
    dimnames = list(paste0("f", 1:10), rownames(locs)))

test_that("neighbour weights follow Banksy's kNN_median kernel", {
    w <- neighborhood_operators(locs, 5L)$w
    expect_equal(unname(Matrix::rowSums(w)), rep(1, n))
    expect_true(all(Matrix::diag(w) == 0))
    expect_equal(unname(Matrix::rowSums(w > 0)), rep(5, n))
    k <- dbscan::kNN(locs, k = 5)
    d <- k$dist[1, ]
    wt <- exp(-d^2 / stats::median(d)^2)
    expect_equal(unname(w[1, k$id[1, ]]), unname(wt / sum(wt)))
})

test_that("reduceData returns the pcaParam list with npcs columns", {
    res <- reduceData(x, neighborhoodPCAParam(lambda = 0.5, k_geom = 5,
        npcs = 4), spatial_locs = locs)
    expect_named(res, c("u", "d", "v", "sdev", "eigenvalues"))
    expect_equal(dim(res$u), c(n, 4L))
    expect_true(any(grepl("_H0$", rownames(res$v))))
})

test_that("lambda = 0 ignores the neighbourhood", {
    res0 <- reduceData(x, neighborhoodPCAParam(lambda = 0, k_geom = 5,
        npcs = 3), spatial_locs = locs)
    expect_false(any(grepl("_H0$", rownames(res0$v))))
})

test_that("matches the Banksy package on expression", {
    skip_if_not_installed("Banksy")
    skip_if_not_installed("SpatialExperiment")
    se <- SpatialExperiment::SpatialExperiment(assays = list(e = x),
        spatialCoords = locs)
    se <- Banksy::computeBanksy(se, assay_name = "e", k_geom = 5,
        verbose = FALSE)
    w <- neighborhood_operators(locs, 5L)$w
    expect_equal(as.matrix(x %*% Matrix::t(w)),
        SummarizedExperiment::assay(se, "H0"), ignore_attr = TRUE)
    se <- Banksy::runBanksyPCA(se, lambda = 0.5, npcs = 4, seed = 1)
    pb <- SingleCellExperiment::reducedDim(se, "PCA_M0_lam0.5")
    pn <- reduceData(x, neighborhoodPCAParam(lambda = 0.5, k_geom = 5,
        npcs = 4), spatial_locs = locs)$u
    expect_equal(abs(diag(stats::cor(pb, pn))), rep(1, 4), tolerance = 1e-6)
})

test_that("runNeighborhoodPCA stores a dimObj and records its param", {
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    g <- runNeighborhoodPCA(g, feats_to_use = "hvf", lambda = 0.8,
        npcs = 5, verbose = FALSE)
    emb <- getDimReduction(g, reduction_method = "neighborhood_pca",
        name = "neighborhood_pca", output = "matrix")
    expect_equal(ncol(emb), 5L)
    expect_setequal(rownames(emb), spatIDs(g))
    h <- g@parameters[[length(g@parameters)]]
    expect_s4_class(attr(h, "attachments")$param, "neighborhoodPCAParam")
})

test_that("the AGF (gradient) block matches the Banksy package", {
    skip_if_not_installed("Banksy")
    skip_if_not_installed("SpatialExperiment")
    se <- SpatialExperiment::SpatialExperiment(assays = list(e = x),
        spatialCoords = locs)
    se <- Banksy::computeBanksy(se, assay_name = "e", k_geom = c(5, 8),
        compute_agf = TRUE, verbose = FALSE)
    h1 <- neighborhood_agf(x, neighborhood_operators(locs, 8L))
    expect_equal(h1, SummarizedExperiment::assay(se, "H1"),
        ignore_attr = TRUE)
    se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = 0.5, npcs = 4,
        seed = 1)
    pb <- SingleCellExperiment::reducedDim(se, "PCA_M1_lam0.5")
    pn <- reduceData(x, neighborhoodPCAParam(lambda = 0.5, k_geom = c(5, 8),
        npcs = 4, compute_agf = TRUE), spatial_locs = locs)$u
    expect_equal(abs(diag(stats::cor(pb, pn))), rep(1, 4), tolerance = 1e-6)
})
