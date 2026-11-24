# runMisty() ####

options(giotto.use_conda = FALSE, giotto.no_python_warn = TRUE)

# mistyR attaches purrr and future when it runs; detach them when this file
# ends so later test files see the usual search path
.misty_search <- search()
withr::defer(for (p in rev(setdiff(search(), .misty_search))) {
    detach(p, character.only = TRUE)
})

# cells in a square with three planted dependences:
#   jx = sum of `a` over the neighbours within 15 (juxtaview)
#   px = gaussian (l = 20) sum of `b` over the other cells (paraview)
#   inn = `cc` of the same cell (intraview)
.misty_sim <- function(n = 400L, side = 220, seed = 11L) {
    GiottoUtils::local_seed(seed)
    xy <- matrix(stats::runif(2L * n, 0, side), n, 2L)
    ids <- sprintf("cell_%03d", seq_len(n))
    d <- as.matrix(stats::dist(xy))
    a <- stats::rnorm(n)
    b <- stats::rnorm(n)
    cc <- stats::rnorm(n)
    dd <- stats::rnorm(n)
    nb <- (d <= 15) * 1
    diag(nb) <- 0
    wp <- exp(-d^2 / 20^2)
    diag(wp) <- 0
    jx <- drop(nb %*% a) + stats::rnorm(n, sd = 0.3)
    px <- drop(wp %*% b) + stats::rnorm(n, sd = 0.3)
    inn <- cc + stats::rnorm(n, sd = 0.2)
    m <- rbind(a = a, b = b, cc = cc, dd = dd, jx = jx, px = px, inn = inn)
    colnames(m) <- ids
    lab <- c("t1", "t2", "t3")[1L + (xy[, 1L] > side / 2) + (xy[, 2L] > side / 2)]
    rna <- matrix(stats::rpois(3L * n, 2), 3L, n,
        dimnames = list(c("g1", "g2", "g3"), ids))
    g <- createGiottoObject(
        expression = list(cell = list(
            prot = list(normalized = m),
            rna = list(normalized = rna))),
        expression_feat = c("prot", "rna"),
        spatial_locs = data.table::data.table(sdimx = xy[, 1L],
            sdimy = xy[, 2L], cell_ID = ids),
        verbose = FALSE)
    addCellMetadata(g, spat_unit = "cell", feat_type = "prot",
        new_metadata = data.table::data.table(cell_ID = ids, label = lab),
        by_column = TRUE, column_cell_ID = "cell_ID")
}

.misty_stat <- function(res, tg, ms) {
    s <- res$improvements.stats
    s$mean[s$target == tg & s$measure == ms]
}

.misty_contrib <- function(res, tg) {
    s <- res$contributions.stats[res$contributions.stats$target == tg, ]
    stats::setNames(s$mean, s$view)
}

test_that("planted juxta, para and intra dependences go to their views", {
    g <- .misty_sim()
    res <- runMisty(g, feat_type = "prot", model = "linear",
        return_gobject = FALSE, verbose = FALSE)
    expect_gt(.misty_stat(res, "jx", "gain.R2"), 30)
    expect_gt(.misty_stat(res, "px", "gain.R2"), 30)
    expect_gt(.misty_stat(res, "inn", "intra.R2"), 80)
    expect_lt(.misty_stat(res, "inn", "gain.R2"), 5)
    expect_lt(.misty_stat(res, "dd", "multi.R2"), 5)
    cj <- .misty_contrib(res, "jx")
    expect_gt(cj[["juxta.15"]], max(abs(cj[c("intra", "para.20")])))
    cp <- .misty_contrib(res, "px")
    expect_gt(cp[["para.20"]], max(abs(cp[c("intra", "juxta.15")])))
    imp <- res$importances.aggregated
    top <- function(v, tg) {
        x <- imp[imp$view == v & imp$Target == tg, ]
        x$Predictor[which.max(x$Importance)]
    }
    expect_identical(top("juxta.15", "jx"), "a")
    expect_identical(top("para.20", "px"), "b")
    expect_identical(top("intra", "inn"), "cc")

    rf <- runMisty(g, feat_type = "prot", targets = c("jx", "px"),
        feats_to_use = c("a", "b", "cc", "dd", "jx", "px"),
        return_gobject = FALSE, verbose = FALSE)
    expect_gt(.misty_stat(rf, "jx", "gain.R2"), 20)
    expect_gt(.misty_stat(rf, "px", "gain.R2"), 20)
})

test_that("results equal mistyR's on the same views", {
    skip_if_not_installed("mistyR")
    skip_if_not_installed("ranger")
    g <- .misty_sim(n = 250L)
    X <- .misty_read_expr(getExpression(g, values = "normalized",
        feat_type = "prot", output = "exprObj")[],
        c("a", "b", "cc", "dd", "jx", "px", "inn"))
    pos <- .misty_positions(g, "cell", NULL, rownames(X), 1, "")
    W <- .misty_weights(g, "cell", pos, 15, NULL, 20, "gaussian", NULL, 0,
        NULL, NULL, 1, FALSE, "")
    views <- .misty_build_views(g, "cell", "prot", X, W, "label", NULL,
        rownames(X), "", FALSE)
    mv <- mistyR::create_initial_view(as.data.frame(views[[1]]$data))
    for (v in views[-1]) {
        mv <- mistyR::add_views(mv, mistyR::create_view(v$name,
            as.data.frame(v$data), v$abbrev))
    }
    key <- function(r) {
        s <- r$improvements.stats
        s <- s[order(s$target, s$measure), ]
        cs <- r$contributions.stats
        cs <- cs[order(cs$target, cs$view), ]
        ia <- r$importances.aggregated
        ia <- ia[order(ia$view, ia$Predictor, ia$Target), ]
        list(s$mean, cs$mean, cs$fraction, cs$p.mean, ia$Importance)
    }
    for (cfg in list(list("linear", FALSE), list("linear", TRUE),
        list("random_forest", FALSE))) {
        dir_m <- tempfile("misty_m")
        dir_g <- tempfile("misty_g")
        on.exit(unlink(c(dir_m, dir_g), recursive = TRUE), add = TRUE)
        suppressWarnings(suppressMessages(if (cfg[[1]] == "linear") {
            mistyR::run_misty(mv, results.folder = dir_m,
                model.function = linear_model, bypass.intra = cfg[[2]])
        } else {
            mistyR::run_misty(mv, results.folder = dir_m,
                model.function = random_forest_model,
                bypass.intra = cfg[[2]])
        }))
        rm <- suppressMessages(mistyR::collect_results(dir_m))
        rg <- runMisty(g, feat_type = "prot", composition = "label",
            model = cfg[[1]], bypass_intra = cfg[[2]],
            results_folder = dir_g, return_gobject = FALSE, verbose = FALSE)
        expect_equal(key(rg), key(rm), tolerance = 1e-6)
        rf <- suppressMessages(mistyR::collect_results(dir_g))
        expect_equal(key(rg), key(rf), tolerance = 1e-10)
        expect_identical(names(rg), names(rm))
        expect_identical(colnames(rg$importances.aggregated),
            colnames(rm$importances.aggregated))
    }
})

test_that("ridge meta-model and p-values match ridge::linearRidge", {
    skip_if_not_installed("ridge")
    GiottoUtils::local_seed(3)
    n <- 300
    X <- matrix(stats::rnorm(n * 4), n, 4)
    X[, 2] <- X[, 1] + stats::rnorm(n, sd = 0.3)
    y <- drop(X %*% c(1, 0.5, -0.3, 0)) + stats::rnorm(n)
    dat <- data.frame(X, y = y)
    fit <- ridge::linearRidge(y ~ ., dat, lambda = "automatic",
        scaling = "corrForm")
    ours <- .misty_ridge(X, y)
    expect_equal(unname(ours$coef), unname(stats::coef(fit)),
        tolerance = 1e-10)
    expect_equal(unname(ours$pval),
        unname(ridge::pvals(fit)$pval[, fit$chosen.nPCs]), tolerance = 1e-8)
})

test_that("folds equal mistyR's caret::createFolds under L'Ecuyer-CMRG", {
    skip_if_not_installed("caret")
    y <- withr::with_seed(1, stats::rnorm(523))
    lecuyer <- function(x, k) {
        withr::with_seed(42, caret::createFolds(x, k = k),
            .rng_kind = "L'Ecuyer-CMRG")
    }
    set.seed(5)
    before <- get(".Random.seed", envir = globalenv())
    expect_identical(.misty_folds(y, 10L, 42L), unname(lecuyer(y, 10)))
    expect_identical(.misty_folds(seq_len(523), 7L, 42L),
        unname(lecuyer(seq_len(523), 7)))
    expect_identical(get(".Random.seed", envir = globalenv()), before)
    expect_identical(RNGkind()[1L], "Mersenne-Twister")
})

test_that("stores R2 in the feature metadata, keeps the caller's RNG", {
    g <- .misty_sim(n = 200L)
    set.seed(99)
    before <- get(".Random.seed", envir = globalenv())
    g2 <- runMisty(g, feat_type = "prot", targets = c("jx", "px"),
        feats_to_use = c("a", "b", "jx", "px"), model = "linear",
        name = "m1", verbose = FALSE)
    expect_identical(get(".Random.seed", envir = globalenv()), before)
    fm <- fDataDT(g2, feat_type = "prot")
    cols <- c("m1_intra_R2", "m1_multi_R2", "m1_gain_R2")
    expect_true(all(cols %in% colnames(fm)))
    expect_true(all(is.na(fm[!feat_ID %in% c("jx", "px"), m1_gain_R2])))
    expect_equal(fm[match(c("jx", "px"), feat_ID), m1_multi_R2] -
        fm[match(c("jx", "px"), feat_ID), m1_intra_R2],
        fm[match(c("jx", "px"), feat_ID), m1_gain_R2])
    g3 <- runMisty(g2, feat_type = "prot", targets = c("jx", "px"),
        feats_to_use = c("a", "b", "jx", "px"), model = "linear",
        name = "m1", verbose = FALSE)
    fm3 <- fDataDT(g3, feat_type = "prot")
    expect_identical(sum(grepl("^m1_", colnames(fm3))), 3L)
    expect_true(any(grepl("_misty", names(g3@parameters))))
})

test_that("cells are matched to positions by ID and views are honoured", {
    g <- .misty_sim(n = 300L)
    args <- list(feat_type = "prot", targets = c("jx", "px"),
        feats_to_use = c("a", "b", "jx", "px"), model = "linear",
        composition = "label", return_gobject = FALSE, verbose = FALSE)
    ref <- do.call(runMisty, c(list(g), args))
    sl <- getSpatialLocations(g, output = "spatLocsObj")
    sl[] <- sl[][rev(seq_len(nrow(sl[])))]
    gs <- setSpatialLocations(g, sl, verbose = FALSE)
    expect_equal(do.call(runMisty, c(list(gs), args)), ref)

    sl_dt <- getSpatialLocations(g, output = "data.table")
    g <- addCellMetadata(g, feat_type = "prot", new_metadata = data.table::data.table(
        cell_ID = sl_dt$cell_ID, left = sl_dt$sdimx < 150), by_column = TRUE,
        column_cell_ID = "cell_ID")
    keep <- sl_dt$cell_ID[sl_dt$sdimx < 150]
    gv <- subset(g, left, view = "left")
    sub <- subsetGiotto(g, cell_ids = keep, verbose = FALSE)
    expect_equal(do.call(runMisty, c(list(gv, view = "left"), args)),
        do.call(runMisty, c(list(sub), args)))
})

test_that("bypass_intra, extra views, nn paraviews and stored networks", {
    g <- .misty_sim(n = 250L)
    base <- list(feat_type = "prot", targets = c("jx", "px"),
        feats_to_use = c("a", "b", "jx", "px"), model = "linear",
        return_gobject = FALSE, verbose = FALSE)
    by <- do.call(runMisty, c(list(g, bypass_intra = TRUE), base))
    ct <- by$contributions[by$contributions$view %in% c("intra", "p.intra"), ]
    expect_identical(ct$value, rep(c(0, 1), 2))
    expect_true(all(by$improvements$value[
        by$improvements$measure == "intra.R2"] == 0))

    ex <- do.call(runMisty, c(list(g, composition = "label",
        extra_views = list(rna = list(feat_type = "rna",
            feats = c("g1", "g2"), views = c("intra", "juxta")))), base))
    expect_setequal(unique(ex$contributions.stats$view),
        c("intra", "juxta.15", "para.20", "juxta.composition.15",
            "para.composition.20", "intra.rna", "juxta.rna.15"))
    expect_setequal(ex$importances.aggregated$Predictor[
        ex$importances.aggregated$view == "juxta.composition.15"],
        c("t1", "t2", "t3"))

    nn <- do.call(runMisty, c(list(g, juxta_radius = NULL, para_nn = 10L),
        base))
    expect_setequal(unique(nn$contributions.stats$view),
        c("intra", "para.20"))

    gn <- createSpatialNetwork(g, method = "Delaunay",
        delaunay_method = "deldir", maximum_distance_delaunay = "auto",
        name = "del", verbose = FALSE)
    ref <- do.call(runMisty, c(list(gn, para_l = NULL), base))
    st <- do.call(runMisty, c(list(gn, para_l = NULL, juxta_network = "del"),
        base))
    expect_equal(st$improvements.stats, ref$improvements.stats)

    two <- do.call(runMisty, c(list(g, para_l = c(10, 30),
        para_family = "exponential"), base))
    expect_setequal(unique(two$contributions.stats$view),
        c("intra", "juxta.15", "para.10", "para.30"))
})

test_that("targets run in parallel on a future plan with the same results", {
    skip_on_cran()
    skip_if_not_installed("future")
    skip_if_not_installed("future.apply")
    g <- .misty_sim(n = 200L)
    args <- list(feat_type = "prot", targets = c("jx", "px", "inn"),
        model = "random_forest", num_trees = 20L, return_gobject = FALSE,
        verbose = FALSE)
    ref <- do.call(runMisty, c(list(g), args))
    old <- future::plan(future::multisession, workers = 2L)
    on.exit(future::plan(old), add = TRUE)
    par <- do.call(runMisty, c(list(g), args))
    expect_equal(par, ref)
})

test_that("input checks", {
    g <- .misty_sim(n = 100L)
    expect_error(runMisty(g, feat_type = "prot", targets = "nope",
        verbose = FALSE), "not features")
    expect_error(runMisty(g, feat_type = "prot", targets = "jx",
        feats_to_use = c("a", "b"), verbose = FALSE), "among")
    expect_error(runMisty(g, feat_type = "prot", composition = "nope",
        model = "linear", verbose = FALSE), "not a cell metadata column")
    expect_error(runMisty(g, feat_type = "prot", coord_scale = 0,
        verbose = FALSE), "> 0")
    expect_error(runMisty(g, feat_type = "prot", para_nn = 5L,
        para_network = "x", verbose = FALSE))
    expect_error(runMisty(g, feat_type = "prot", model = "linear",
        extra_views = list(r = list(feat_type = "rna", feats = "zz")),
        verbose = FALSE), "not features")
})

test_that("sample_fraction keeps the default forest and accepts a smaller share", {
    skip_if_not_installed("ranger")
    g <- .misty_sim(n = 200L)
    args <- list(feat_type = "prot", targets = c("jx", "px"),
        feats_to_use = c("a", "b", "jx", "px"), num_trees = 20L, cores = 1L,
        return_gobject = FALSE, verbose = FALSE)
    a <- do.call(runMisty, c(list(g), args))
    expect_identical(do.call(runMisty, c(list(g, sample_fraction = 1), args)), a)
    s <- do.call(runMisty, c(list(g, sample_fraction = 0.3), args))
    expect_false(identical(s$improvements, a$improvements))
    expect_true(all(is.finite(s$improvements$value)))
    expect_error(do.call(runMisty, c(list(g, sample_fraction = 2), args)))
})

.misty_add_fov <- function(g, side = 220) {
    sl <- getSpatialLocations(g, output = "data.table")
    addCellMetadata(g, feat_type = "prot", new_metadata = data.table::data.table(
        cell_ID = sl$cell_ID,
        fov = 1L + (sl$sdimx > side / 2) + 2L * (sl$sdimy > side / 2)),
        by_column = TRUE, column_cell_ID = "cell_ID")
}

test_that("sample_by equals mistyR run per sample on whole-slide views", {
    skip_if_not_installed("mistyR")
    g <- .misty_add_fov(.misty_sim(n = 400L))
    X <- .misty_read_expr(getExpression(g, values = "normalized",
        feat_type = "prot", output = "exprObj")[],
        c("a", "b", "cc", "dd", "jx", "px", "inn"))
    pos <- .misty_positions(g, "cell", NULL, rownames(X), 1, "")
    W <- .misty_weights(g, "cell", pos, 15, NULL, 20, "gaussian", NULL, 0,
        NULL, NULL, 1, FALSE, "")
    views <- .misty_build_views(g, "cell", "prot", X, W, NULL, NULL,
        rownames(X), "", FALSE)
    fov <- pDataDT(g, feat_type = "prot")[match(rownames(X), cell_ID), fov]
    dirs <- file.path(withr::local_tempdir(), paste0("fov", 1:4))
    for (s in 1:4) {
        r <- which(fov == s)
        mv <- mistyR::create_initial_view(as.data.frame(
            views[[1]]$data[r, , drop = FALSE]))
        for (v in views[-1]) {
            mv <- mistyR::add_views(mv, mistyR::create_view(v$name,
                as.data.frame(v$data[r, , drop = FALSE]), v$abbrev))
        }
        suppressWarnings(suppressMessages(mistyR::run_misty(mv,
            results.folder = dirs[s], model.function = linear_model)))
    }
    rm <- suppressMessages(mistyR::collect_results(dirs))
    rg <- runMisty(g, feat_type = "prot", sample_by = "fov", min_cells = 50L,
        model = "linear", return_gobject = FALSE, verbose = FALSE)
    key <- function(r) {
        s <- r$improvements.stats
        s <- s[order(s$target, s$measure), ]
        cs <- r$contributions.stats
        cs <- cs[order(cs$target, cs$view), ]
        ia <- r$importances.aggregated
        ia <- ia[order(ia$view, ia$Predictor, ia$Target), ]
        list(s$mean, s$sd, cs$mean, cs$p.mean, ia$Importance, ia$nsamples)
    }
    expect_equal(key(rg), key(rm), tolerance = 1e-6)
    expect_setequal(unique(rg$improvements$sample), as.character(1:4))
    expect_identical(nrow(rg$improvements), nrow(rm$improvements))
})

test_that("sample_by leaves out small samples and constant targets", {
    g <- .misty_add_fov(.misty_sim(n = 300L))
    cm <- pDataDT(g, feat_type = "prot")
    lab <- ifelse(cm$fov == 4L, "tiny", as.character(cm$fov))
    lab[cm$fov == 4L][-(1:5)] <- "3"
    g <- addCellMetadata(g, feat_type = "prot",
        new_metadata = data.table::data.table(cell_ID = cm$cell_ID,
            smp = lab), by_column = TRUE, column_cell_ID = "cell_ID")
    ex <- getExpression(g, values = "normalized", feat_type = "prot")
    m <- ex[]
    m["dd", cm$cell_ID[cm$fov == 1L]] <- 1
    ex[] <- m
    g <- setExpression(g, ex, verbose = FALSE)
    args <- list(feat_type = "prot", targets = c("jx", "dd"),
        feats_to_use = c("a", "b", "jx", "dd"), model = "linear",
        sample_by = "smp", min_cells = 20L, return_gobject = FALSE)
    expect_message(res <- do.call(runMisty, c(list(g), args)),
        "1 with fewer than 20 cells")
    expect_setequal(unique(res$improvements$sample), c("1", "2", "3"))
    fitted <- unique(res$improvements[, c("sample", "target")])
    expect_false(any(fitted$sample == "1" & fitted$target == "dd"))
    expect_identical(nrow(fitted), 5L)
    expect_identical(res$importances.aggregated$nsamples[
        res$importances.aggregated$Target == "dd"][1L], 2L)
    expect_error(do.call(runMisty, c(list(g), utils::modifyList(args,
        list(min_cells = 1000L, verbose = FALSE)))), "no .smp. group")
})

test_that("sample_by runs samples in parallel and writes one folder each", {
    skip_on_cran()
    skip_if_not_installed("future")
    skip_if_not_installed("future.apply")
    skip_if_not_installed("mistyR")
    g <- .misty_add_fov(.misty_sim(n = 300L))
    args <- list(feat_type = "prot", targets = c("jx", "px"),
        feats_to_use = c("a", "b", "jx", "px"), num_trees = 20L,
        sample_by = "fov", min_cells = 30L, return_gobject = FALSE,
        verbose = FALSE)
    out <- withr::local_tempdir()
    ref <- do.call(runMisty, c(list(g, results_folder = out), args))
    expect_setequal(list.files(out), as.character(1:4))
    files <- suppressMessages(mistyR::collect_results(file.path(out, 1:4)))
    expect_equal(files$improvements.stats$mean[order(
        files$improvements.stats$target, files$improvements.stats$measure)],
        ref$improvements.stats$mean[order(ref$improvements.stats$target,
            ref$improvements.stats$measure)], tolerance = 1e-10)
    seq <- do.call(runMisty, c(list(g), args))
    old <- future::plan(future::multisession, workers = 2L)
    on.exit(future::plan(old), add = TRUE)
    expect_equal(do.call(runMisty, c(list(g), args)), seq)
})
