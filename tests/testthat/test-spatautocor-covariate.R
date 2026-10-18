# The covariate null: Moran's I of within-group residuals and the conditioned
# join count. Checks: one group reduces to the plain null, the residual cross
# term matches a dense computation, and on a library-size gradient the
# covariate null stops flagging genes that only follow library size.

.cv_sim <- function(N = 3000L, seed = 5L) {
    set.seed(seed)
    xy <- cbind(runif(N, 0, 30), runif(N, 0, 30))
    nn <- dbscan::kNN(xy, k = 6)$id
    a <- rep(seq_len(N), 6); b <- as.vector(nn)
    e <- unique(data.frame(from = pmin(a, b), to = pmax(a, b)))
    ids <- paste0("c", seq_len(N))
    s <- exp(log(3) * xy[, 1] / 30 + rnorm(N, 0, 0.2))
    s <- s / mean(s)
    blob <- as.numeric(sqrt((xy[, 1] - 15)^2 + (xy[, 2] - 15)^2) < 5)
    # spatial genes are a small share, so the blob barely moves library size
    n_lib <- 120L; n_spat <- 12L
    mu <- c(rep(c(0.3, 1, 2), length.out = n_lib),
            rep(c(0.5, 1, 2), length.out = n_spat))
    rate <- cbind(outer(s, mu[seq_len(n_lib)]),
                  outer(s * (1 + 3 * blob), mu[n_lib + seq_len(n_spat)]))
    counts <- t(matrix(rpois(length(rate), rate), N))
    lib <- colSums(counts)
    x <- log1p(t(t(counts) * (1e3 / pmax(lib, 1))))
    dimnames(x) <- list(c(paste0("lib_", seq_len(n_lib)),
                          paste0("spat_", seq_len(n_spat))), ids)
    list(x = Matrix::Matrix(x, sparse = TRUE),
         net = data.frame(from = ids[e$from], to = ids[e$to]),
         lib = stats::setNames(lib, ids))
}
.cv_bins <- function(v, K = 10L) {
    stats::setNames(cut(v, unique(quantile(v, 0:K / K)),
                        include.lowest = TRUE, labels = FALSE), names(v))
}

test_that("one group reduces the covariate null to the plain one", {
    skip_if_not_installed("dbscan")
    f <- .cv_sim(N = 600L)
    one <- stats::setNames(rep("all", ncol(f$x)), colnames(f$x))
    st <- GiottoClass::analyzeData(f$x, GiottoClass::edgeStatsParam(),
        network = f$net, groups = one)
    expect_equal(.moran_covariate_test(st)[, -"feats"], .moran_test(st)[, -"feats"])

    b <- processData(f$x, thresholdParam("binarize", threshold = "kmeans"))
    sb <- GiottoClass::analyzeData(b, GiottoClass::edgeStatsParam(),
        network = f$net, groups = one)
    cov <- .joincount_covariate_test(sb)
    perm <- .joincount_test(sb)
    expect_equal(cov$expected, perm$expected, tolerance = 0.01)
    expect_equal(cov$variance, perm$variance, tolerance = 0.02)
})

test_that("the residual Moran cross term matches a dense computation", {
    skip_if_not_installed("dbscan")
    f <- .cv_sim(N = 400L)
    grp <- .cv_bins(f$lib, 4L)
    st <- GiottoClass::analyzeData(f$x, GiottoClass::edgeStatsParam(),
        network = f$net, groups = grp)
    got <- .moran_covariate_test(st)
    X <- as.matrix(f$x)
    R <- X - t(apply(X, 1, function(v) ave(v, grp)))
    N <- ncol(X)
    W <- matrix(0, N, N)
    W[cbind(match(f$net$from, colnames(X)), match(f$net$to, colnames(X)))] <- 1
    W <- W + t(W)
    ref <- N / sum(W) * rowSums((R %*% W) * R) / rowSums(R^2)
    expect_equal(got$statistic, unname(ref))
})

test_that("the covariate null stops flagging library-size-only genes", {
    skip_if_not_installed("dbscan")
    f <- .cv_sim()
    grp <- .cv_bins(f$lib)
    flag <- function(res) {
        sig <- stats::p.adjust(res$p.value, "BH") < 0.05
        c(lib = mean(sig[startsWith(res$feats, "lib_")]),
          spat = mean(sig[startsWith(res$feats, "spat_")]))
    }
    auroc <- function(res) {
        pos <- startsWith(res$feats, "spat_")
        r <- rank(res$z)
        (sum(r[pos]) - sum(pos) * (sum(pos) + 1) / 2) / (sum(pos) * sum(!pos))
    }
    # moran keeps its power; binSpect's binarized test pays some at this size
    min_power <- c(moran = 0.5, binspect = 0.25)
    for (method in c("moran", "binspect")) {
        perm <- analyzeData(f$x, spatAutoCorParam(method), network = f$net)
        cov <- analyzeData(f$x, spatAutoCorParam(method, null = "covariate"),
            network = f$net, groups = grp)
        lab <- function(x) paste(method, x)
        expect_gt(flag(perm)[["lib"]], 0.3, label = lab("permutation"))
        expect_lt(flag(cov)[["lib"]], 0.1, label = lab("covariate"))
        expect_gt(flag(cov)[["spat"]], min_power[[method]], label = lab("power"))
    }
    cov <- analyzeData(f$x, spatAutoCorParam("moran", null = "covariate"),
        network = f$net, groups = grp)
    expect_gt(auroc(cov), 0.95)
})

test_that("the giotto method derives library-size groups", {
    skip_if_not_installed("GiottoData")
    g <- suppressMessages(suppressWarnings(
        GiottoData::loadGiottoMini("visium", verbose = FALSE)))
    res <- analyzeData(g, spatAutoCorParam("moran", null = "covariate"))
    expect_true(all(is.finite(res$z)))
    expect_error(spatAutoCorParam("geary", null = "covariate"), "moran")
    expect_error(analyzeData(g, spatAutoCorParam("moran", null = "covariate",
        covariate = "nope")), "metadata")
})
