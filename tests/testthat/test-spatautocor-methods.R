# The spatAutoCorParam compositions: each method equals its blocks run by
# hand, the weight-matrix path equals terra's statistic, and the giotto
# entry points run end to end on the visium mini object.

.sm_fixture <- function(G = 5L, N = 70L, seed = 3L) {
    set.seed(seed)
    xy <- cbind(runif(N), runif(N))
    nn <- dbscan::kNN(xy, k = 4)$id
    a <- rep(seq_len(N), 4); b <- as.vector(nn)
    e <- unique(data.frame(from = pmin(a, b), to = pmax(a, b)))
    ids <- paste0("c", seq_len(N))
    lam <- cbind(8 * xy[, 1], matrix(3, N, G - 1))
    m <- Matrix::Matrix(t(matrix(rpois(N * G, lam), N)), sparse = TRUE) * 1
    dimnames(m) <- list(paste0("g", seq_len(G)), ids)
    list(m = m, net = data.frame(from = ids[e$from], to = ids[e$to]), xy = xy)
}

test_that("moran and geary params are the edge stats plus the statistic", {
    skip_if_not_installed("dbscan")
    f <- .sm_fixture()
    st <- GiottoClass::analyzeData(f$m, GiottoClass::edgeStatsParam(),
                                   network = f$net)
    for (case in list(list("moran", .moran_test), list("geary", .geary_test))) {
        got <- analyzeData(f$m, spatAutoCorParam(case[[1]]), network = f$net)
        ref <- case[[2]](st)
        ref <- ref[match(got$feats, ref$feats)]
        expect_equal(got$statistic, ref$statistic)
        expect_equal(got$adj.p.value,
                     p.adjust(ref$p.value, "fdr")[seq_len(nrow(ref))])
        expect_false(is.unsorted(rev(got$z)))
    }
    expect_equal(analyzeData(f$m, spatAutoCorParam("moran"), network = f$net)$feats[1],
                 "g1")
})

test_that("the binspect param binarizes, then counts joins", {
    skip_if_not_installed("dbscan")
    f <- .sm_fixture()
    for (rule in c("kmeans", "rank")) {
        b <- processData(f$m, thresholdParam("binarize", threshold = rule))
        ref <- .joincount_test(GiottoClass::analyzeData(b,
            GiottoClass::edgeStatsParam(), network = f$net))
        got <- analyzeData(f$m, spatAutoCorParam("binspect", bin_method = rule),
                           network = f$net)
        ref <- ref[match(got$feats, ref$feats)]
        expect_equal(got$statistic, ref$statistic)
        expect_equal(got$n_on, ref$n_on)
    }
})

test_that("the weight-matrix path equals terra's statistic for any W", {
    skip_if_not_installed("dbscan")
    f <- .sm_fixture()
    N <- ncol(f$m)
    W <- matrix(0, N, N)
    W[cbind(match(f$net$from, colnames(f$m)), match(f$net$to, colnames(f$m)))] <- 1
    W <- W + t(W)
    W <- W / rowSums(W)  # row-standardised, so asymmetric
    vals <- t(as.matrix(f$m))
    for (method in c("moran", "geary")) {
        got <- .spat_autocor_weights(vals, colnames(vals), W, method)
        ref <- apply(vals, 2, function(v) terra::autocor(v, W, method))
        expect_equal(got$statistic[match(names(ref), got$feats)], unname(ref))
    }
})

test_that("binSpect fast and the giotto method agree on a mini object", {
    skip_if_not_installed("GiottoData")
    g <- suppressMessages(suppressWarnings(
        GiottoData::loadGiottoMini("visium", verbose = FALSE)))
    p <- spatAutoCorParam("binspect")
    res <- analyzeData(g, p, spatial_network_name = "Delaunay_network")
    expect_equal(nrow(res), nrow(getExpression(g, values = "normalized",
        output = "matrix")))
    expect_equal(binSpect(g, implementation = "fast", verbose = FALSE), res)

    g2 <- binSpect(g, implementation = "fast", return_gobject = TRUE)
    fm <- getFeatureMetadata(g2, output = "data.table")
    expect_true(all(c("binSpect.pval", "binspect.z", "binspect.padj") %in% names(fm)))
    expect_equal(fm$binSpect.pval[match(res$feats, fm$feat_ID)], res$adj.p.value)
})

test_that("spatialAutoCor is the analyzeData call, and joincount is binspect", {
    skip_if_not_installed("GiottoData")
    g <- suppressMessages(suppressWarnings(
        GiottoData::loadGiottoMini("visium", verbose = FALSE)))
    expect_equal(spatialAutoCor(g, method = "moran", return_gobject = FALSE),
                 analyzeData(g, spatAutoCorParam("moran")))
    expect_equal(spatialAutoCor(g, method = "joincount", return_gobject = FALSE),
                 spatialAutoCor(g, method = "binspect", return_gobject = FALSE))
    expect_s4_class(spatAutoCorParam("joincount"), "binSpectParam")
    g2 <- spatialAutoCor(g, method = "joincount")
    fm <- getFeatureMetadata(g2, output = "data.table")
    expect_true(all(c("binspect.z", "binspect.padj", "binSpect.pval") %in% names(fm)))
})
