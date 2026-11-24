# The HMRF engine and its smfishHmrf equivalence are tested in spatialHMRF;
# these tests cover Giotto's verb, wrapper and composition.
options(giotto.has_conda = FALSE, giotto.use_conda = FALSE,
    giotto.no_python_warn = TRUE)

sim_tissue <- function(n = 300, seed = 1) {
    set.seed(seed)
    xy <- matrix(runif(2 * n), ncol = 2) * 10
    dom <- 1L + (xy[, 1] > 5)
    y <- matrix(rnorm(n * 6), n) + c(-1.5, 1.5)[dom]
    ids <- paste0("c", seq_len(n))
    rownames(y) <- ids
    knn <- dbscan::kNN(xy, k = 5)$id
    e <- unique(t(apply(cbind(rep(seq_len(n), 5), as.vector(knn)), 1, sort)))
    list(y = y, dom = dom,
        net = data.frame(from = ids[e[, 1]], to = ids[e[, 2]]))
}

test_that("clusterData(hmrfParam) returns one domain per cell", {
    skip_if_not_installed("spatialHMRF")
    s <- sim_tissue()
    res <- clusterData(s$y, hmrfParam(k = 2, beta = 2, nstart = 10,
        zscore = "none"), network = s$net)
    expect_named(res, c("cell_ID", "cluster"))
    expect_equal(res$cell_ID, rownames(s$y))
    expect_equal(dim(attr(res, "prob")), c(nrow(s$y), 2L))
    expect_true(is.list(attr(res, "blocks")))
    expect_gt(aricode::ARI(res$cluster, s$dom), 0.9)
})

test_that("clusterData(hmrfParam) equals spatialHMRF::hmrf()", {
    skip_if_not_installed("spatialHMRF")
    s <- sim_tissue()
    res <- clusterData(s$y, hmrfParam(k = 2, beta = 2, nstart = 10,
        zscore = "none"), network = s$net)
    fit <- spatialHMRF::hmrf(s$y, s$net, k = 2, beta = 2, zscore = "none",
        nstart = 10)
    expect_identical(res$cluster, fit$class)
    expect_equal(attr(res, "prob"), fit$prob)
})

test_that("hmrfParam passes the start and stopping options to spatialHMRF", {
    skip_if_not_installed("spatialHMRF")
    s <- sim_tissue()
    p <- hmrfParam(k = 2, beta = 2, nstart = 5, zscore = "none",
        subsample = 100, label_tol = 0)
    res <- clusterData(s$y, p, network = s$net)
    fit <- spatialHMRF::hmrf(s$y, s$net, k = 2, beta = 2, zscore = "none",
        nstart = 5, subsample = 100, label_tol = 0)
    expect_identical(res$cluster, fit$class)
    expect_error(hmrfParam(k = 4, subsample = 2))
})

test_that("runHMRF stores domains, probabilities and its param", {
    skip_if_not_installed("spatialHMRF")
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    g <- runHMRF(g, feats_to_use = "hvf", k = 4, betas = c(0, 2),
        nstart = 5, verbose = FALSE)
    md <- pDataDT(g)
    expect_true(all(c("hmrf_k4_b0", "hmrf_k4_b2") %in% names(md)))
    enr <- getSpatialEnrichment(g, name = "hmrf_k4_b2", output = "data.table")
    expect_equal(nrow(enr), nrow(md))
    h <- g@parameters[[length(g@parameters)]]
    expect_s4_class(attr(h, "attachments")$param, "hmrfParam")
})

test_that("runHMRF fits a k x beta grid and hmrfSummary compares it", {
    skip_if_not_installed("spatialHMRF")
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    g <- runHMRF(g, feats_to_use = "hvf", k = c(3, 4), betas = c(0, 4),
        nstart = 5, verbose = FALSE)
    s <- hmrfSummary(g)
    expect_equal(nrow(s), 4L)
    expect_equal(s$k, c(3, 3, 4, 4))
    expect_true(all(s$confidence > 0 & s$confidence <= 1))
    expect_gt(s[k == 3 & beta == 4]$coherence, s[k == 3 & beta == 0]$coherence)
})

test_that("findSpatialDomains passes wrapper arguments through", {
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    g <- createSpatialNetwork(g, method = "kNN", k = 15, name = "knn15",
        verbose = FALSE)
    g <- findSpatialDomains(g, method = "neighborhood",
        spatial_network_name = "knn15", npcs = 10, name = "dom",
        verbose = FALSE)
    expect_false(anyNA(pDataDT(g)$dom))
    h <- g@parameters[[length(g@parameters)]]
    expect_equal(attr(h, "attachments")$param$npcs, 10L)
})

test_that("runHMRF and findSpatialDomains fit a spatial enrichment", {
    skip_if_not_installed("spatialHMRF", minimum_version = "0.1.1")
    skip_if_not_installed("GiottoData")
    g <- GiottoData::loadGiottoMini("visium")
    sv <- spatialAutoCor(g, method = "moran", return_gobject = FALSE)
    g <- spatialCoexpression(g, feats = head(sv$feats, 200), n_modules = 8,
        verbose = FALSE)
    g <- runHMRF(g, spat_enr_name = "spatial_module", k = 4, nstart = 5,
        verbose = FALSE)
    expect_false(anyNA(pDataDT(g)$hmrf_k4_b2))
    h <- g@parameters[[length(g@parameters)]]
    expect_equal(attr(h, "attachments")$param@param$zscore, "feature")
    g <- findSpatialDomains(g, method = "hmrf",
        spat_enr_name = "spatial_module", k = 4, nstart = 5, name = "dom")
    expect_equal(sort(unique(pDataDT(g)$dom)), 1:4)
    h <- g@parameters[[length(g@parameters)]]
    expect_equal(attr(h, "attachments")$param@param$zscore, "feature")
    expect_equal(attr(h, "attachments")$param@param$nstart, 5L)
})
