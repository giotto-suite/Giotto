# Spatial autocorrelation statistics from the edgeStatsParam contract. Ground
# truth: a dense-matrix Moran's I, and the mean and variance over real
# permutations of the values across cells.

.sac_fixture <- function(G = 4L, N = 80L, seed = 11L) {
    set.seed(seed)
    xy <- cbind(runif(N), runif(N))
    nn <- dbscan::kNN(xy, k = 4)$id
    a <- rep(seq_len(N), 4); b <- as.vector(nn)
    e <- unique(data.frame(from = pmin(a, b), to = pmax(a, b)))
    ids <- paste0("c", seq_len(N))
    # feature 1 follows x, the rest are noise; counts with zeros
    lam <- cbind(10 * xy[, 1], matrix(3, N, G - 1))
    m <- Matrix::Matrix(t(matrix(rpois(N * G, lam), N)), sparse = TRUE) * 1
    dimnames(m) <- list(paste0("g", seq_len(G)), ids)
    net <- data.frame(from = ids[e$from], to = ids[e$to],
                      weight = round(runif(nrow(e), 0.5, 2), 2))
    list(m = m, net = net, e = e)
}
.sac_stats <- function(m, net, weighted = FALSE) {
    GiottoClass::analyzeData(m, GiottoClass::edgeStatsParam(weighted = weighted),
                             network = net)
}
.perm_moments <- function(f, test, nperm = 1500L, weighted = FALSE,
                          binarize = FALSE) {
    m <- f$m
    if (binarize) m <- (m > 3) * 1
    sims <- replicate(nperm, {
        p <- sample(ncol(m))
        mp <- m[, p]
        colnames(mp) <- colnames(m)
        test(.sac_stats(mp, f$net, weighted))$statistic
    })
    an <- test(.sac_stats(m, f$net, weighted))
    list(an = an, mean = rowMeans(sims), sd = apply(sims, 1, sd))
}

test_that("Moran's I equals the dense formula, weighted or not", {
    skip_if_not_installed("dbscan")
    f <- .sac_fixture()
    for (wt in c(FALSE, TRUE)) {
        w <- if (wt) f$net$weight else rep(1, nrow(f$net))
        W <- matrix(0, ncol(f$m), ncol(f$m))
        W[cbind(f$e$from, f$e$to)] <- w
        W <- W + t(W)
        z <- as.matrix(f$m) - Matrix::rowMeans(f$m)
        ref <- ncol(f$m) / sum(W) * rowSums((z %*% W) * z) / rowSums(z^2)
        got <- .moran_test(.sac_stats(f$m, f$net, wt))
        expect_equal(got$statistic, unname(ref))
    }
    expect_lt(got$p.value[1], 1e-6)  # the gradient feature
})

test_that("Moran, Geary and join-count moments match real permutations", {
    skip_if_not_installed("dbscan")
    set.seed(2)
    f <- .sac_fixture()
    for (case in list(list(.moran_test, FALSE, FALSE),
                      list(.geary_test, TRUE, FALSE),
                      list(.joincount_test, FALSE, TRUE))) {
        r <- .perm_moments(f, case[[1]], weighted = case[[2]],
                           binarize = case[[3]])
        expect_equal(r$an$expected, r$mean, tolerance = 0.02)
        expect_equal(sqrt(r$an$variance), r$sd, tolerance = 0.08)
    }
})

test_that("the join count is the on-on edge count with its 2x2 odds ratio", {
    skip_if_not_installed("dbscan")
    f <- .sac_fixture()
    b <- (f$m > 3) * 1
    got <- .joincount_test(.sac_stats(b, f$net))
    bd <- as.matrix(b)
    on_a <- bd[, f$e$from, drop = FALSE]; on_b <- bd[, f$e$to, drop = FALSE]
    n11 <- rowSums(on_a * on_b); n00 <- rowSums((1 - on_a) * (1 - on_b))
    n10 <- rowSums(on_a * (1 - on_b)); n01 <- rowSums((1 - on_a) * on_b)
    expect_equal(got$statistic, unname(n11))
    expect_equal(got$n_on, unname(rowSums(bd)))
    expect_equal(got$odds_ratio, unname(n11 * n00 / ((n10 + n01) / 2)^2))
})

test_that("ACAT combines p-values and keeps a shared p", {
    expect_equal(.acat(matrix(0.2, 3, 4)), rep(0.2, 3))
    expect_lt(.acat(matrix(c(1e-8, 0.5, 0.5), 1)), 1e-7)
})

test_that("spatAutoCorParam builds the concrete class per method", {
    expect_s4_class(spatAutoCorParam("moran"), "moranParam")
    expect_s4_class(spatAutoCorParam("geary"), "gearyParam")
    p <- spatAutoCorParam("binspect", bin_method = "rank", percentage_rank = 20)
    expect_s4_class(p, "binSpectParam")
    expect_true(is(p, "spatAutoCorParam"))
    expect_equal(p$percentage_rank, 20)
    expect_null(spatAutoCorParam("moran")$bin_method)
})
