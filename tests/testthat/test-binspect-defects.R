# Fixes to the legacy binSpect implementations.

.bd_fixture <- function(N = 40L, G = 6L, seed = 2L) {
    set.seed(seed)
    ids <- paste0("c", seq_len(N))
    expr <- matrix(rpois(G * N, 3), G, N, dimnames = list(paste0("g", 1:G), ids))
    bin <- (expr > 3) * 1
    a <- rep(seq_len(N), each = 3L)
    b <- ((a - 1L + rep(c(1L, 2L, 5L), N)) %% N) + 1L
    e <- unique(data.table::data.table(i = pmin(a, b), j = pmax(a, b)))
    net <- data.table::data.table(from = ids[e$i], to = ids[e$j],
                                  distance = runif(nrow(e), 1, 2))
    net$distance[1] <- 0  # two coincident cells
    sym <- rbind(net, data.table::data.table(from = net$to, to = net$from,
                                             distance = net$distance))
    list(expr = expr, bin = bin, net = sym)
}

test_that("the matrix implementation runs serially and keeps zero-distance edges", {
    f <- .bd_fixture()
    run <- function(impl, ...) {
        r <- binSpectSingleMatrix(f$expr, spatial_network = f$net,
            bin_matrix = f$bin, implementation = impl, get_av_expr = FALSE,
            get_high_expr = FALSE, verbose = FALSE, ...)
        r[order(feats)]
    }
    dt <- run("data.table", cores = 1)
    mx <- run("matrix", do_parallel = FALSE)
    expect_equal(mx$p.value, dt$p.value)
    expect_equal(mx$estimate, dt$estimate)
})

test_that("binSpectMulti combines p-values with 2k degrees of freedom", {
    f <- .bd_fixture()
    nets <- list(f$net, f$net[-(1:6)])
    res <- binSpectMultiMatrix(f$expr, spatial_networks = nets,
        bin_method = "rank", get_av_expr = FALSE, get_high_expr = FALSE,
        verbose = FALSE)
    # rank: one trial per network and percentage_rank (default c(10, 30))
    comb <- res$combined[, list(x = -2 * sum(log(p.value)), k = .N), by = feats]
    expect_equal(unique(comb$k), length(nets) * 2L)
    ref <- stats::pchisq(comb$x, df = 2 * comb$k, lower.tail = FALSE)
    expect_equal(res$simple$p.val[match(comb$feats, res$simple$feats)], ref)
})

test_that("binSpect with spatial_network_k can return the gobject", {
    skip_if_not_installed("GiottoData")
    g <- suppressMessages(suppressWarnings(
        GiottoData::loadGiottoMini("visium", verbose = FALSE)))
    g2 <- binSpect(g, spatial_network_k = c(5, 8), bin_method = "rank",
                   return_gobject = TRUE, verbose = FALSE)
    fm <- getFeatureMetadata(g2, output = "data.table")
    expect_true("binSpect.pval" %in% names(fm))
    expect_false(all(is.na(fm$binSpect.pval)))
})

test_that("a constant feature is called off instead of failing kmeans", {
    m <- rbind(a = rep(2, 30), b = c(rep(0, 15), rep(5, 15)))
    for (algo in c("kmeans", "kmeans_arma", "kmeans_arma_subset")) {
        b <- kmeans_binarize_wrapper(m, kmeans_algo = algo, seed = 1,
                                     extreme_nr = 5, sample_nr = 5)
        expect_equal(unname(b[1, ]), rep(0, 30), label = algo)
        expect_equal(unname(b[2, ]), rep(c(0, 1), each = 15), label = algo)
    }
})

test_that("the kmeans_arma subset is seeded and works on few cells", {
    set.seed(9)
    x <- c(rpois(200, 1), rpois(100, 8))
    s1 <- .kmeans_arma_subset_binarize(x, extreme_nr = 20, sample_nr = 50, seed = 3)
    s2 <- .kmeans_arma_subset_binarize(x, extreme_nr = 20, sample_nr = 50, seed = 3)
    expect_identical(s1, s2)
    expect_length(.kmeans_arma_subset_binarize(x[1:40], extreme_nr = 20,
        sample_nr = 50, seed = 3), 40L)
})
