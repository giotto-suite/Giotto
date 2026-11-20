options(giotto.use_conda = FALSE, giotto.no_python_warn = TRUE,
    giottodisk.par_workers = 2L, gkernels.n_threads = 2L)


make_synth <- function(n_per_fov = 400L, seed = 1L) {
    set.seed(seed)
    fovs <- rep(1:3, each = n_per_fov)
    n <- length(fovs)
    fov_bg <- c(1, 2, 1.5)[fovs]
    cell_bg <- exp(stats::rnorm(n, 0, 0.2))
    bgl <- 6 * fov_bg * cell_bg
    ln <- function(s) exp(stats::rnorm(n, 0, s))
    pos <- function(p, level) ifelse(stats::runif(n) < p, level * ln(0.3), 0)
    iso <- cbind(`Ms IgG1` = bgl * ln(0.15), `Rb IgG` = 1.1 * bgl * ln(0.15))
    tg <- cbind(
        A = 2 * bgl * ln(0.2) + pos(0.3, 400),
        B = 3 * bgl * ln(0.2) + pos(0.2, 800),
        C = 1.5 * bgl * ln(0.2) + pos(0.4, 300),
        D = 1.2 * bgl * ln(0.2) + pos(0.1, 600),
        Bg1 = 1.1 * bgl * ln(0.2),
        Bg2 = 1.3 * bgl * ln(0.2))
    ids <- sprintf("c_1_%d_%d", fovs, seq_len(n))
    tg <- round(tg, 2); iso <- round(iso, 2)
    rownames(tg) <- rownames(iso) <- ids
    list(tg = tg, iso = iso, fov = fovs, ids = ids)
}

build_g <- function(s, backend = NULL) {
    ex <- list(protein = list(raw = t(s$tg)), protein_neg = list(raw = t(s$iso)))
    sl <- data.table::data.table(cell_ID = s$ids, sdimx = seq_along(s$ids), sdimy = s$fov * 100)
    args <- list(expression = ex, spatial_locs = sl, verbose = FALSE)
    if (!is.null(backend)) args$backend <- backend
    do.call(GiottoClass::createGiottoObject, args)
}

norm_mat <- function(g, name) {
    t(as.matrix(.read_expr_matrix(GiottoClass::getExpression(g,
        feat_type = "protein", values = name, output = "matrix"))))
}

norm <- function(g, norm_methods, name = NULL, cofactor = 5,
    clip_quantile = c(0.01, 0.99), expression_values = "raw") {
    if (is.null(name)) {
        return(normalizeGiotto(g, feat_type = "protein", expression_values = expression_values,
            norm_methods = norm_methods, cofactor = cofactor, clip_quantile = clip_quantile,
            verbose = FALSE))
    }
    normalizeGiotto(g, feat_type = "protein", expression_values = expression_values,
        norm_methods = norm_methods, cofactor = cofactor, clip_quantile = clip_quantile,
        name = name, verbose = FALSE)
}

S <- make_synth()
G <- build_g(S)
X <- t(S$tg)

test_that("gamma_pearson equals per-feature Gamma GLM Pearson residuals", {
    got <- t(.gamma_pearson_norm(X, c(0.01, 0.99)))
    lcs <- pmax(log(rowSums(S$tg)), 0)
    ref <- sapply(colnames(S$tg), function(f) {
        x <- S$tg[, f]
        y <- pmin(pmax(x, stats::quantile(x, 0.01)), stats::quantile(x, 0.99))
        stats::residuals(stats::glm(y ~ offset(lcs), family = stats::Gamma(link = "log"),
            control = stats::glm.control(epsilon = 1e-14, maxit = 100)), type = "pearson")
    })
    expect_equal(unname(got), unname(ref), tolerance = 1e-8)
})

test_that("clr is the centred log-ratio across the features of each cell", {
    ref <- t(apply(S$tg, 1, function(x) log1p(x / exp(sum(log1p(x)) / length(x)))))
    expect_equal(t(.clr_norm(X)), ref, tolerance = 1e-12)
    expect_equal(processData(X, normParam("clr")), .clr_norm(X))
    expect_equal(as.matrix(processData(methods::as(X, "dgCMatrix"), normParam("clr"))), .clr_norm(X),
        tolerance = 1e-12)
})

test_that("normalizeGiotto stores each method under 'normalized' or `name`", {
    for (m in c("clr", "gamma_pearson", "arcsinh")) {
        y <- norm_mat(norm(G, norm_methods = m), "normalized")
        expect_identical(dim(y), dim(S$tg), info = m)
        expect_setequal(rownames(y), S$ids)
        expect_setequal(colnames(y), colnames(S$tg))
        expect_false(anyNA(y))
        expect_true(all(is.finite(norm_mat(norm(G, norm_methods = m, name = m), m))), info = m)
    }
    expect_equal(norm_mat(norm(G, norm_methods = "arcsinh"), "normalized")[S$ids, colnames(S$tg)],
        asinh(S$tg / 5))
    expect_equal(norm_mat(norm(G, norm_methods = "clr"), "normalized")[S$ids, colnames(S$tg)],
        t(.clr_norm(X)), tolerance = 1e-12)
})

test_that("the arcsinh cofactor takes one value or a value per feature", {
    f <- colnames(S$tg)
    one <- stats::setNames(150, f[1])
    y <- norm_mat(norm(G, norm_methods = "arcsinh", cofactor = c(5, one)), "normalized")
    expect_equal(y[S$ids, f], asinh(sweep(S$tg, 2L, c(150, rep(5, length(f) - 1L)), "/")))
    all_named <- stats::setNames(seq_along(f) + 1, f)
    expect_equal(processData(X, normParam("arcsinh", c = all_named)), asinh(X / all_named))
    expect_equal(processData(X, normParam("arcsinh", c = 10)), asinh(X / 10))
    expect_error(norm(G, norm_methods = "arcsinh", cofactor = one), "no value for")
    expect_error(norm(G, norm_methods = "arcsinh", cofactor = c(5, nope = 2)), "not features: nope")
    expect_error(norm(G, norm_methods = "arcsinh", cofactor = c(5, 6)), "one cofactor")
    expect_error(norm(G, norm_methods = "arcsinh", cofactor = c(5, 6, A = 2)), "more than one unnamed")
    expect_error(norm(G, norm_methods = "arcsinh", cofactor = 0), "> 0")
})

test_that("log1p is normParam('log') with base e", {
    expect_equal(processData(X, normParam("log", base = exp(1))), log1p(X), tolerance = 1e-12)
})

test_that("gamma_pearson clips at the quantiles, or not with c(0, 1)", {
    s <- exp(pmax(log(colSums(X)), 0))
    rel <- t(t(X) / s)
    expect_equal(.gamma_pearson_norm(X, c(0, 1)), rel / rowMeans(rel) - 1, tolerance = 1e-12)
    xc <- pmin(X, apply(X, 1, stats::quantile, 0.9, names = FALSE))
    rel <- t(t(xc) / s)
    expect_equal(.gamma_pearson_norm(X, 0.9), rel / rowMeans(rel) - 1, tolerance = 1e-12)
    y <- norm_mat(norm(G, norm_methods = "gamma_pearson", clip_quantile = c(0, 1)), "normalized")
    expect_equal(y[S$ids, colnames(S$tg)], t(.gamma_pearson_norm(X, c(0, 1))), tolerance = 1e-12)
    expect_error(norm(G, norm_methods = "gamma_pearson", clip_quantile = c(0.9, 0.1)), "lower, upper")
})

test_that("a non-default spat_unit is normalized like the default one", {
    gn <- GiottoClass::createGiottoObject(expression = list(nucleus = list(
        protein = list(raw = t(S$tg)))), verbose = FALSE)
    gn <- normalizeGiotto(gn, spat_unit = "nucleus", feat_type = "protein",
        norm_methods = "clr", verbose = FALSE)
    yn <- t(as.matrix(GiottoClass::getExpression(gn, spat_unit = "nucleus", feat_type = "protein",
        values = "normalized", output = "matrix")))
    expect_equal(yn[S$ids, ], norm_mat(norm(G, norm_methods = "clr"), "normalized")[S$ids, ])
})

test_that("non-raw input and name = input are errors", {
    g <- norm(G, norm_methods = "gamma_pearson", name = "gp")
    expect_error(norm(g, norm_methods = "arcsinh", expression_values = "gp"), "raw intensities")
    s2 <- S
    s2$tg[1L, "A"] <- Inf
    expect_error(norm(build_g(s2), norm_methods = "clr"), "finite")
    expect_error(norm(G, norm_methods = "clr", name = "raw"), "overwritten")
})

test_that("a feature that is 0 in every cell gives 0, not non-finite values", {
    s2 <- S
    s2$tg[, "D"] <- 0
    expect_warning(g <- norm(build_g(s2), norm_methods = "gamma_pearson"), "0 in every cell are set to 0: D")
    y <- norm_mat(g, "normalized")
    expect_true(all(is.finite(y)))
    expect_true(all(y[, "D"] == 0))
})

test_that("a GiottoDisk-backed object gives the in-memory result and writes a store", {
    skip_if_not_installed("GiottoDisk")
    vault <- tempfile("pn_vault_")
    on.exit(unlink(vault, recursive = TRUE), add = TRUE)
    gd <- build_g(S, backend = vault)
    expect_true(inherits(GiottoClass::getExpression(gd, feat_type = "protein", values = "raw",
        output = "matrix"), "parquetExprBase"))
    for (m in c("gamma_pearson", "clr", "arcsinh")) {
        gd <- norm(gd, norm_methods = m, name = m)
        st <- GiottoClass::getExpression(gd, feat_type = "protein", values = m, output = "matrix")
        expect_true(inherits(st, "parquetExprBase"), info = m)
        mem <- norm_mat(norm(G, norm_methods = m, name = m), m)[S$ids, colnames(S$tg)]
        expect_equal(norm_mat(gd, m)[S$ids, colnames(S$tg)], mem, tolerance = 1e-12, info = m)
    }
    keep <- S$ids[S$fov != 2][seq(1, 800, by = 3)]
    gs <- norm(GiottoClass::subsetGiotto(gd, cell_ids = keep, verbose = FALSE),
        norm_methods = "gamma_pearson", name = "sub")
    gm <- norm(GiottoClass::subsetGiotto(G, cell_ids = keep, verbose = FALSE),
        norm_methods = "gamma_pearson", name = "sub")
    ys <- norm_mat(gs, "sub")
    expect_identical(nrow(ys), length(keep))
    expect_equal(ys[keep, colnames(S$tg)], norm_mat(gm, "sub")[keep, colnames(S$tg)], tolerance = 1e-12)
})
