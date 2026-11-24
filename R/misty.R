# MISTy multiview spatial models ####

#' @title Multiview spatial models of feature expression (MISTy)
#' @name runMisty
#' @description
#' Fits MISTy, the multiview model of Tanevski et al. (2022, *Genome
#' Biology* 23:97), on the cells of a Giotto object. For every target feature
#' it measures how much of the target's variation is explained by the other
#' features of the same cell (intraview), by the sum over the cell's
#' immediate neighbours (juxtaview), by a distance-weighted sum over its wider
#' neighbourhood (paraview) and, optionally, by the cell-type composition of
#' both neighbourhoods and by other feature types or enrichment scores. It is
#' the algorithm of the mistyR package (`run_misty()` with
#' `random_forest_model` or `linear_model`, then `collect_results()`): the
#' same view models, seeds, cross-validation, meta-model, p-values and result
#' tables, so mistyR's `plot_*()` functions draw its output. One model is
#' fitted on the whole slide, or one per sample (e.g. per FOV) with
#' `sample_by`.
#'
#' It runs on all cells, also of whole slides and GiottoDisk-backed objects:
#' only the features used are read, cells are matched to their positions by
#' cell ID, neighbourhoods come from Giotto spatial networks and every view is
#' a sparse product (cost linear in the number of cells and edges), and the
#' meta-model's p-values are computed in closed form (mistyR's
#' `ridge::pvals()` forms an N x N matrix, which limits it to about 20,000
#' cells).
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type of the targets and of the intra-, juxta- and
#'   paraview
#' @param expression_values expression values of `feat_type` to model
#' @param targets features of `feat_type` to model. `NULL` (default): all of
#'   `feats_to_use`, or all features when both are `NULL`.
#' @param feats_to_use features of `feat_type` in the intraview, juxtaview and
#'   paraview, i.e. the predictors. Must include `targets`. `NULL` (default):
#'   `targets`.
#' @param spat_loc_name name of the spatial locations. `NULL`: the default.
#' @param coord_scale numeric > 0. Factor that converts the spatial
#'   coordinates to the units of `juxta_radius`, `para_l`, `para_cutoff` and
#'   `para_zoi` (e.g. 0.12028 for CosMx pixels to micrometres).
#' @param juxta_radius numeric or `NULL`. Juxtaview: sum over the Delaunay
#'   neighbours (`deldir`, as mistyR) no farther than this distance. `NULL`:
#'   no juxtaview.
#' @param juxta_network name of a stored spatial network of `spat_unit` to use
#'   for the juxtaview instead of a Delaunay network; its edges no longer than
#'   `juxta_radius` are kept.
#' @param para_l numeric (one or more) or `NULL`. Paraview: sum over the other
#'   cells weighted by `para_family` with length scale `para_l`, one paraview
#'   per value. `NULL`: no paraview.
#' @param para_family weight of a cell at distance d: `"gaussian"`
#'   `exp(-d^2 / l^2)`, `"exponential"` `exp(-d / l)`, `"linear"`
#'   `max(0, 1 - d / l)` or `"constant"` 1 (with `para_nn` nearest
#'   neighbours, `round(l)` when `para_nn` is `NULL`), as mistyR's
#'   `add_paraview()`.
#' @param para_cutoff numeric or `NULL`. Cells farther than
#'   `para_cutoff * para_l` are left out of the paraview (mistyR sums over all
#'   cells). `NULL`: 4 for `"gaussian"` and 16 for `"exponential"` (weights
#'   below 1.2e-7 are dropped), 1 for `"linear"` (exact).
#' @param para_zoi numeric. Zone of indifference: cells closer than this get
#'   weight 0 in the paraview (mistyR's `zoi`).
#' @param para_nn integer or `NULL`. Limit the paraview to the `para_nn`
#'   nearest neighbours of each cell (mistyR's `nn`).
#' @param para_network name of a stored spatial network of `spat_unit` to use
#'   for the paraview neighbours instead of building one. Weights come from the
#'   distances; a kNN network is used as directed.
#' @param composition `NULL` or a cell metadata column with a label per cell
#'   (e.g. cell types). Adds the juxtaview and paraview(s) of the one-hot
#'   labels: how many cells of each label are next to, or around, the cell.
#' @param extra_views `NULL` or a named list of further predictor sources.
#'   Each element is a list with either `feat_type` (+ optional
#'   `expression_values`, default `"normalized"`, and `feats`, default all;
#'   at most 500 features) or `enrichment` (the name of a spatial enrichment,
#'   + optional `feat_type`), and `views`, a subset of
#'   `c("intra", "juxta", "para")` (default `c("juxta", "para")`): the source's
#'   values in the cell itself, summed over the juxta neighbours, or weighted
#'   over the paraview.
#' @param sample_by `NULL` (default) or a cell metadata column (e.g. `"fov"`).
#'   `NULL` fits one model per target on all cells (the whole slide); a column
#'   fits one model per target in each of its groups (samples), on views built
#'   from the whole slide. See Details.
#' @param min_cells integer. With `sample_by`, samples with fewer cells (or
#'   fewer than `cv_folds`) are left out.
#' @param model view model: `"random_forest"` (ranger, mistyR's
#'   `random_forest_model`: 100 trees, impurity importance, out-of-bag
#'   predictions) or `"linear"` (mistyR's `linear_model`: least squares,
#'   out-of-fold predictions, coefficients as importances).
#' @param bypass_intra logical. Leave the intraview out of the meta-model
#'   (mistyR's `bypass.intra`): the intraview model becomes a constant, so
#'   multi.R2 is what the other views explain on their own.
#' @param cv_folds integer. Folds of the cross-validated performance estimate
#'   and of the linear view models' out-of-fold predictions.
#' @param num_trees,mtry random forest: number of trees and of features tried
#'   at each split (`NULL`: ranger's default, the square root of the number of
#'   predictors).
#' @param sample_fraction random forest: share of the cells drawn (with
#'   replacement) for each tree. `NULL` (default) draws as many cells as there
#'   are, as ranger and mistyR do; a smaller share such as 0.1 fits large data
#'   several times faster.
#' @param cores integer. Threads of each random forest. `NA`: Giotto's
#'   `giotto.cores` option. With a parallel `future` plan the targets (or,
#'   with `sample_by`, the samples) are instead fitted in parallel, one thread
#'   each; every worker loads the packages, so on one machine the threads
#'   usually match the workers' speed with less memory (per FOV on 146,000
#'   cells: 7.8 min and 6.2 GB with 4 threads, 8.5 min and 11.7 GB with 4
#'   workers).
#' @param seed integer. Seed of the random forests and the folds (mistyR's
#'   default 42). The caller's random number state is restored on exit.
#' @param results_folder `NULL` or a folder to write mistyR's result files
#'   to (`coefficients.txt`, `performance.txt`, `importances_*.txt`), which
#'   `mistyR::collect_results()` reads. With `sample_by`, one subfolder per
#'   sample.
#' @param name prefix of the feature metadata columns
#' @param view character. Name of a `GiottoClass::giottoView` slotted on
#'   `gobject`. Supplied, the object is narrowed once before the expression,
#'   spatial locations, cell metadata and enrichments are read, so the model
#'   uses the view's cells.
#' @param return_gobject logical. Return the giotto object with the results in
#'   the feature metadata (default) or the result list.
#' @param verbose be verbose
#' @returns With `return_gobject = TRUE`, the giotto object with
#' `<name>_intra_R2`, `<name>_multi_R2` and `<name>_gain_R2` (in %, the mean
#' over the cross-validation folds, and over the samples with `sample_by`)
#' added to the feature metadata of `feat_type` (`NA` for features that are
#' not targets). Otherwise a list shaped as `mistyR::collect_results()`
#' returns it (data.frames `improvements`, `contributions`, `importances`,
#' `improvements.stats`, `contributions.stats`, `importances.aggregated`),
#' which mistyR's `plot_improvement_stats()`, `plot_view_contributions()`,
#' `plot_interaction_heatmap()` and `plot_interaction_communities()` take.
#' Its `sample` column is the value of `sample_by`; for a whole-slide model,
#' the absolute path of `results_folder`, or `name`.
#' @details
#' For each target and each view a model predicts the target from the view's
#' features: out-of-bag predictions for random forests, out-of-fold
#' predictions for linear models. Like mistyR, a view built from
#' `feats_to_use` (intra-, juxta- and paraview) leaves out the target's own
#' column; composition and extra views keep all their columns. The meta-model
#' regresses the target on the views' predictions: ordinary least squares
#' when it has one predictor, otherwise ridge regression with the penalty
#' chosen automatically (Cule and De Iorio 2012, as `ridge::linearRidge()`
#' with `lambda = "automatic"`, which mistyR uses). Its coefficients are the
#' view contributions, with p-values of a normal test (ridge, as
#' `ridge::pvals()`) or a t test (least squares). Performance is
#' estimated by `cv_folds`-fold cross-validation of the meta-model (folds
#' stratified by the target's quantiles, as `caret::createFolds()`):
#' intra.R2 is the R2 of the intraview alone, multi.R2 that of all views,
#' gain.R2 = multi.R2 - intra.R2 (in %; negative fold values are set to 0, as
#' mistyR does), with RMSE versions and one-sided t tests of the gain
#' (`p.R2`, `p.RMSE`). Importances are standardised per target and view and
#' weighted by 1 - p-value of the view's contribution, as `collect_results()`
#' does.
#'
#' Whole slide or per sample: with `sample_by`, the views are built once from
#' the whole slide (a cell's neighbours across a sample border still count)
#' and the models are fitted on each sample's cells, as mistyR's
#' `run_misty()` on each sample followed by `collect_results()`. The
#' `*.stats` tables then hold means and standard deviations over the
#' samples, and `improvements` one row per sample, so the gains can be mapped
#' back onto the tissue. Per-sample models are MISTy's own design: each learns
#' the relationships within one region, so differences between regions do not
#' count as a neighbourhood effect, but each sees fewer cells (lower
#' intra.R2, noisier for rare targets). A target that is constant in a sample
#' is not fitted there.
#'
#' Differences from mistyR: views are built from Giotto networks (Delaunay by
#' `deldir`, `dbscan` radius or kNN searches) on positions matched by cell ID;
#' the paraview is truncated at `para_cutoff * para_l`; the linear view models
#' use `cv_folds` folds (mistyR's `linear_model()` always uses 10); the
#' result files of a run with only the intraview have one column fewer than
#' mistyR's. With the same views and settings the R2, contributions and
#' importances equal mistyR's.
#'
#' Cost: views are sparse products, linear in the number of cells and edges
#' (about 130 edges per cell for a gaussian paraview with `para_l` = 20 um at
#' cell densities of tissue). The models hold the views as dense matrices
#' (cells x features of all views); a random forest with 100 trees on 146,000
#' cells and 21 predictors takes about 11 s with 4 threads, so a run costs
#' about (targets x views) of those. Extra views are capped at 500 features
#' each: select features (e.g. highly variable genes) or use enrichment
#' scores.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' feats <- head(fDataDT(g)[order(-perc_cells), feat_ID], 10)
#'
#' res <- runMisty(g, targets = feats[1:2], feats_to_use = feats,
#'     juxta_radius = 300, para_l = 300, model = "linear",
#'     return_gobject = FALSE)
#' res$improvements.stats
#' @md
#' @export
runMisty <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    targets = NULL,
    feats_to_use = NULL,
    spat_loc_name = NULL,
    coord_scale = 1,
    juxta_radius = 15,
    juxta_network = NULL,
    para_l = 20,
    para_family = c("gaussian", "exponential", "linear", "constant"),
    para_cutoff = NULL,
    para_zoi = 0,
    para_nn = NULL,
    para_network = NULL,
    composition = NULL,
    extra_views = NULL,
    sample_by = NULL,
    min_cells = 200L,
    model = c("random_forest", "linear"),
    bypass_intra = FALSE,
    cv_folds = 10L,
    num_trees = 100L,
    mtry = NULL,
    sample_fraction = NULL,
    cores = NA,
    seed = 42L,
    results_folder = NULL,
    name = "misty",
    view = NULL,
    return_gobject = TRUE,
    verbose = NULL) {
    fname <- "[runMisty]"
    checkmate::assert_class(gobject, "giotto")
    model <- match.arg(model)
    para_family <- match.arg(para_family)
    checkmate::assert_character(targets, min.len = 1L, any.missing = FALSE,
        unique = TRUE, null.ok = TRUE)
    checkmate::assert_character(feats_to_use, min.len = 1L,
        any.missing = FALSE, unique = TRUE, null.ok = TRUE)
    checkmate::assert_string(spat_loc_name, null.ok = TRUE)
    checkmate::assert_number(coord_scale, lower = 0, finite = TRUE)
    checkmate::assert_number(juxta_radius, lower = 0, finite = TRUE,
        null.ok = TRUE)
    checkmate::assert_string(juxta_network, null.ok = TRUE)
    checkmate::assert_numeric(para_l, lower = 0, finite = TRUE,
        any.missing = FALSE, min.len = 1L, unique = TRUE, null.ok = TRUE)
    checkmate::assert_number(para_cutoff, lower = 0, finite = TRUE,
        null.ok = TRUE)
    checkmate::assert_number(para_zoi, lower = 0, finite = TRUE)
    checkmate::assert_count(para_nn, positive = TRUE, null.ok = TRUE)
    checkmate::assert_string(para_network, null.ok = TRUE)
    checkmate::assert_string(composition, null.ok = TRUE)
    checkmate::assert_list(extra_views, names = "unique", null.ok = TRUE)
    checkmate::assert_string(sample_by, null.ok = TRUE)
    checkmate::assert_count(min_cells, positive = TRUE)
    checkmate::assert_flag(bypass_intra)
    checkmate::assert_int(cv_folds, lower = 2L)
    checkmate::assert_count(num_trees, positive = TRUE)
    checkmate::assert_count(mtry, positive = TRUE, null.ok = TRUE)
    checkmate::assert_number(sample_fraction, lower = 1e-6, upper = 1,
        null.ok = TRUE)
    checkmate::assert_int(seed)
    checkmate::assert_string(results_folder, null.ok = TRUE)
    checkmate::assert_string(name, min.chars = 1L)
    checkmate::assert_flag(return_gobject)
    if (coord_scale <= 0) {
        stop(fname, " `coord_scale` must be > 0", call. = FALSE)
    }
    if (any(para_l <= 0)) {
        stop(fname, " `para_l` must be > 0", call. = FALSE)
    }
    package_check("caret", repository = "CRAN")
    if (model == "random_forest") {
        package_check("ranger", repository = "CRAN")
    }
    cores <- GiottoUtils::determine_cores(cores)
    GiottoUtils::local_seed(seed)

    spat_unit <- GiottoClass::set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- GiottoClass::set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )
    values <- match.arg(expression_values,
        unique(c("normalized", "scaled", "custom", expression_values)))
    gobj <- gobject
    if (!is.null(view)) {
        checkmate::assert_string(view)
        gobj <- GiottoClass::resolveRecipe(gobj, view = view,
            slots = c("expression", "spatial_locs", "cell_metadata",
                "spatial_enrichment"),
            spat_unit = spat_unit, feat_type = feat_type)
    }

    ex <- GiottoClass::getExpression(gobj, values = values,
        spat_unit = spat_unit, feat_type = feat_type, output = "exprObj")
    all_feats <- rownames(ex[])
    if (is.null(targets)) targets <- feats_to_use %null% all_feats
    if (is.null(feats_to_use)) feats_to_use <- targets
    miss <- setdiff(c(feats_to_use, targets), all_feats)
    if (length(miss)) {
        stop(GiottoUtils::wrap_txt(fname, "not features of", sQuote(feat_type),
            ":", paste(utils::head(miss, 5L), collapse = ", ")), call. = FALSE)
    }
    if (!all(targets %in% feats_to_use)) {
        stop(fname, " `targets` must be among `feats_to_use`", call. = FALSE)
    }
    X <- .misty_read_expr(ex[], feats_to_use)
    cells <- rownames(X)
    n <- length(cells)
    if (is.null(sample_by)) {
        if (n < cv_folds) {
            stop(fname, " fewer cells (", n, ") than `cv_folds`",
                call. = FALSE)
        }
        .misty_check_targets(X[, targets, drop = FALSE], cv_folds, fname)
        samples <- NULL
    } else {
        if (any(!is.finite(X[, targets]))) {
            stop(fname, " the targets have missing or infinite values",
                call. = FALSE)
        }
        samples <- .misty_samples(gobj, spat_unit, feat_type, sample_by,
            cells, max(min_cells, cv_folds), verbose, fname)
    }

    pos <- .misty_positions(gobj, spat_unit, spat_loc_name, cells,
        coord_scale, fname)
    W <- .misty_weights(gobj, spat_unit, pos, juxta_radius, juxta_network,
        para_l, para_family, para_cutoff, para_zoi, para_nn, para_network,
        coord_scale, verbose, fname)

    views <- .misty_build_views(gobj, spat_unit, feat_type, X, W, composition,
        extra_views, cells, fname, verbose)
    rm(W)
    abbrev <- vapply(views, `[[`, "", "abbrev")
    if (!is.null(results_folder)) {
        .misty_check_names(c(targets, abbrev, names(samples)), fname)
    }
    GiottoUtils::vmsg(.v = verbose, sprintf(
        "%s %d targets, %d cells, views: %s", fname, length(targets), n,
        paste(sprintf("%s (%d)", abbrev,
            vapply(views, function(v) ncol(v$data), 1L)), collapse = ", ")))

    t0 <- Sys.time()
    if (is.null(sample_by)) {
        fits <- list(.misty_fit_sample(X, views, targets, model,
            bypass_intra, cv_folds, num_trees, mtry, sample_fraction, seed,
            cores, parallel = TRUE, verbose = verbose, fname = fname))
        names(fits) <- name
    } else {
        fits <- .misty_fit_samples(samples, X, views, targets, model,
            bypass_intra, cv_folds, num_trees, mtry, sample_fraction, seed,
            cores, verbose, fname)
    }
    rm(views)
    GiottoUtils::vmsg(.v = verbose, sprintf("%s models fitted in %.1f s",
        fname, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    .misty_fit_notes(fits, targets, verbose, fname)

    if (!is.null(results_folder)) {
        dirs <- if (is.null(sample_by)) {
            results_folder
        } else {
            file.path(results_folder, names(fits))
        }
        for (s in which(lengths(fits) > 0L)) {
            .misty_write(dirs[s], names(fits[[s]]), abbrev, fits[[s]])
        }
        if (is.null(sample_by)) names(fits) <- normalizePath(results_folder)
        GiottoUtils::vmsg(.v = verbose, fname, "results written to",
            normalizePath(results_folder))
    }
    res <- .misty_collect(fits, abbrev)

    if (!isTRUE(return_gobject)) return(res)

    st <- data.table::as.data.table(res$improvements.stats)
    st <- data.table::dcast(st[st$measure %in%
        c("intra.R2", "multi.R2", "gain.R2")], target ~ measure,
        value.var = "mean")
    cols <- paste0(name, c("_intra_R2", "_multi_R2", "_gain_R2"))
    new_meta <- data.table::data.table(feat_ID = st$target)
    new_meta[, (cols) := list(st$intra.R2, st$multi.R2, st$gain.R2)]
    gobject <- .misty_set_feat_meta(gobject, spat_unit, feat_type, new_meta)
    GiottoClass::update_giotto_params(gobject, description = "_misty")
}


# internals ####

.misty_read_expr <- function(m, feats) {
    m <- t(as.matrix(.read_expr_matrix(m, feats)))
    storage.mode(m) <- "double"
    m[, feats, drop = FALSE]
}

.misty_check_targets <- function(Y, cv_folds, fname) {
    if (any(!is.finite(Y))) {
        stop(fname, " the targets have missing or infinite values",
            call. = FALSE)
    }
    sds <- apply(Y, 2L, stats::sd)
    if (any(sds == 0)) {
        stop(GiottoUtils::wrap_txt(fname, "targets with zero variance:",
            paste(colnames(Y)[sds == 0], collapse = ", ")), call. = FALSE)
    }
    nu <- apply(Y, 2L, function(v) length(unique(v)))
    if (any(nu < cv_folds)) {
        warning(GiottoUtils::wrap_txt(fname, "targets",
            paste(colnames(Y)[nu < cv_folds], collapse = ", "),
            "have fewer unique values than cv_folds"), call. = FALSE)
    }
    invisible(TRUE)
}

.misty_samples <- function(gobject, spat_unit, feat_type, sample_by, cells,
    min_cells, verbose, fname) {
    lab <- .misty_cell_column(gobject, spat_unit, feat_type, sample_by, cells,
        fname)
    samples <- split(seq_along(cells), factor(lab, levels = sort(unique(lab))))
    small <- lengths(samples) < min_cells
    if (all(small)) {
        stop(GiottoUtils::wrap_txt(fname, "no", sQuote(sample_by),
            "group has", min_cells, "cells or more"), call. = FALSE)
    }
    GiottoUtils::vmsg(.v = verbose, sprintf(paste(
        "%s %d samples by %s (median %d cells); %d with fewer than %d cells",
        "left out"), fname, sum(!small), sample_by,
        as.integer(stats::median(lengths(samples[!small]))), sum(small),
        min_cells))
    samples[!small]
}

.misty_check_names <- function(x, fname) {
    bad <- x[grepl("[[:space:]/\\\\]", x)]
    if (length(bad)) {
        stop(GiottoUtils::wrap_txt(fname, "names with spaces or slashes",
            "cannot be written to mistyR's result files:",
            paste(bad, collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}

.misty_positions <- function(gobject, spat_unit, spat_loc_name, cells,
    coord_scale, fname) {
    sl <- GiottoClass::getSpatialLocations(gobject, spat_unit = spat_unit,
        name = spat_loc_name, output = "data.table")
    i <- match(cells, sl$cell_ID)
    if (anyNA(i)) {
        stop(GiottoUtils::wrap_txt(fname, sum(is.na(i)),
            "cells of the expression have no spatial location, e.g.",
            paste(utils::head(cells[is.na(i)], 3L), collapse = ", ")),
            call. = FALSE)
    }
    pos <- cbind(sl$sdimx[i], sl$sdimy[i]) * coord_scale
    dimnames(pos) <- list(cells, c("sdimx", "sdimy"))
    pos
}

.misty_stored_edges <- function(gobject, spat_unit, net_name, cells, fname) {
    sn <- GiottoClass::getSpatialNetwork(gobject, spat_unit = spat_unit,
        name = net_name, output = "spatialNetworkObj")
    dt <- GiottoClass::getSpatialNetwork(gobject, spat_unit = spat_unit,
        name = net_name, output = "networkDT")
    if (!all(c("from", "to") %in% colnames(dt))) {
        stop(fname, " spatial network ", sQuote(net_name),
            " has no from/to columns", call. = FALSE)
    }
    i <- match(as.character(dt$from), cells)
    j <- match(as.character(dt$to), cells)
    keep <- !is.na(i) & !is.na(j) & i != j
    list(i = i[keep], j = j[keep],
        directed = identical(tolower(methods::slot(sn, "method")), "knn"))
}

.misty_edges_from_dt <- function(dt, cells) {
    i <- match(as.character(dt$from), cells)
    j <- match(as.character(dt$to), cells)
    keep <- !is.na(i) & !is.na(j) & i != j
    list(i = i[keep], j = j[keep])
}

.misty_dist <- function(pos, i, j) {
    sqrt((pos[i, 1L] - pos[j, 1L])^2 + (pos[i, 2L] - pos[j, 2L])^2)
}

# one row per unordered pair, both directions in the matrix
.misty_sym <- function(i, j, w, n) {
    a <- pmin(i, j)
    b <- pmax(i, j)
    u <- !duplicated(a * (n + 1) + b)
    Matrix::sparseMatrix(i = c(a[u], b[u]), j = c(b[u], a[u]),
        x = c(w[u], w[u]), dims = c(n, n))
}

.misty_family_weight <- function(d, family, l, zoi) {
    w <- switch(family,
        gaussian = exp(-d^2 / l^2),
        exponential = exp(-d / l),
        linear = pmax(0, 1 - d / l),
        constant = rep(1, length(d))
    )
    w[d < zoi] <- 0
    w
}

.misty_weights <- function(gobject, spat_unit, pos, juxta_radius,
    juxta_network, para_l, para_family, para_cutoff, para_zoi, para_nn,
    para_network, coord_scale, verbose, fname) {
    n <- nrow(pos)
    cells <- rownames(pos)
    out <- list(juxta = NULL, juxta_radius = juxta_radius, para = list())
    if (!is.null(juxta_radius)) {
        e <- if (is.null(juxta_network)) {
            dt <- GiottoClass::createNetwork(pos,
                GiottoClass::delaunayNetworkParam(method = "deldir",
                    maximum_distance = juxta_radius, include_weight = FALSE,
                    output = "data.table"))
            .misty_edges_from_dt(dt, cells)
        } else {
            .misty_stored_edges(gobject, spat_unit, juxta_network, cells,
                fname)
        }
        d <- .misty_dist(pos, e$i, e$j)
        k <- d <= juxta_radius
        out$juxta <- .misty_sym(e$i[k], e$j[k], rep(1, sum(k)), n)
        GiottoUtils::vmsg(.v = verbose, sprintf(
            paste("%s juxtaview: %.2f neighbours per cell within %g,",
                "%d cells without"),
            fname, Matrix::nnzero(out$juxta) / n, juxta_radius,
            sum(Matrix::rowSums(out$juxta) == 0)))
    }
    if (is.null(para_l)) return(out)

    cut_mult <- para_cutoff %null% c(gaussian = 4, exponential = 16,
        linear = 1, constant = Inf)[[para_family]]
    nn_l <- if (!is.null(para_nn)) {
        rep(as.integer(para_nn), length(para_l))
    } else if (para_family == "constant") {
        as.integer(round(para_l))
    }
    rk <- NULL
    if (!is.null(para_network)) {
        if (!is.null(para_nn)) {
            stop(fname, " give `para_nn` or `para_network`, not both",
                call. = FALSE)
        }
        e <- .misty_stored_edges(gobject, spat_unit, para_network, cells,
            fname)
        d <- .misty_dist(pos, e$i, e$j)
    } else if (!is.null(nn_l)) {
        e <- .misty_edges_from_dt(GiottoClass::createNetwork(pos,
            GiottoClass::kNNNetworkParam(k = min(max(nn_l), n - 1L),
                include_weight = FALSE, output = "data.table")), cells)
        e$directed <- TRUE
        d <- .misty_dist(pos, e$i, e$j)
        o <- order(e$i, d)
        rk <- integer(length(d))
        rk[o] <- sequence(rle(e$i[o])$lengths)
    } else {
        e <- .misty_edges_from_dt(GiottoClass::createNetwork(pos,
            GiottoClass::radiusNetworkParam(eps = cut_mult * max(para_l),
                include_weight = FALSE, output = "data.table")), cells)
        e$directed <- FALSE
        d <- .misty_dist(pos, e$i, e$j)
    }
    for (k in seq_along(para_l)) {
        l <- para_l[k]
        keep <- d <= cut_mult * l
        if (!is.null(rk)) keep <- keep & rk <= nn_l[k]
        w <- .misty_family_weight(d[keep], para_family, l, para_zoi)
        W <- if (isTRUE(e$directed)) {
            Matrix::sparseMatrix(i = e$i[keep], j = e$j[keep], x = w,
                dims = c(n, n))
        } else {
            .misty_sym(e$i[keep], e$j[keep], w, n)
        }
        out$para[[as.character(l)]] <- W
        GiottoUtils::vmsg(.v = verbose, sprintf(
            "%s paraview l = %g: %.1f neighbours per cell", fname, l,
            Matrix::nnzero(W) / n))
    }
    out
}

.misty_onehot <- function(lab) {
    f <- factor(lab)
    m <- Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
        dims = c(length(f), nlevels(f)))
    colnames(m) <- levels(f)
    m
}

.misty_cell_column <- function(gobject, spat_unit, feat_type, column, cells,
    fname) {
    cm <- GiottoClass::list_cell_metadata(gobject, spat_unit = spat_unit)
    fts <- unique(c(feat_type, cm$feat_type))
    fts <- fts[fts %in% cm$feat_type]
    for (ft in fts) {
        md <- GiottoClass::getCellMetadata(gobject, spat_unit = spat_unit,
            feat_type = ft, output = "data.table")
        if (column %in% colnames(md)) {
            v <- md[[column]][match(cells, md$cell_ID)]
            if (anyNA(v)) {
                stop(GiottoUtils::wrap_txt(fname, sQuote(column),
                    "is missing for", sum(is.na(v)), "cells"), call. = FALSE)
            }
            return(v)
        }
    }
    stop(fname, " ", sQuote(column), " is not a cell metadata column",
        call. = FALSE)
}

.misty_extra_source <- function(gobject, spat_unit, feat_type, src, nm,
    cells, fname) {
    checkmate::assert_list(src, names = "unique", .var.name = nm)
    checkmate::assert_subset(names(src),
        c("feat_type", "expression_values", "feats", "enrichment", "views"),
        .var.name = nm)
    vw <- src$views %null% c("juxta", "para")
    checkmate::assert_subset(vw, c("intra", "juxta", "para"), empty.ok = FALSE,
        .var.name = paste0(nm, "$views"))
    if (!is.null(src$enrichment)) {
        checkmate::assert_string(src$enrichment)
        en <- GiottoClass::getSpatialEnrichment(gobject, spat_unit = spat_unit,
            feat_type = src$feat_type %null% feat_type,
            name = src$enrichment, output = "data.table")
        i <- match(cells, en$cell_ID)
        if (anyNA(i)) {
            stop(fname, " enrichment ", sQuote(src$enrichment), " misses ",
                sum(is.na(i)), " cells", call. = FALSE)
        }
        m <- as.matrix(en[i, setdiff(colnames(en), "cell_ID"), with = FALSE])
        rownames(m) <- cells
    } else {
        checkmate::assert_string(src$feat_type, .var.name = paste0(nm,
            "$feat_type"))
        vals <- src$expression_values %null% "normalized"
        ex <- GiottoClass::getExpression(gobject, values = vals,
            spat_unit = spat_unit, feat_type = src$feat_type,
            output = "exprObj")[]
        feats <- src$feats %null% rownames(ex)
        miss <- setdiff(feats, rownames(ex))
        if (length(miss)) {
            stop(GiottoUtils::wrap_txt(fname, nm, ": not features of",
                sQuote(src$feat_type), ":",
                paste(utils::head(miss, 5L), collapse = ", ")), call. = FALSE)
        }
        if (length(feats) > 500L) {
            stop(GiottoUtils::wrap_txt(fname, nm, "has", length(feats),
                "features; extra views take at most 500 (the models hold",
                "every view as a dense cells x features matrix). Select",
                "features, e.g. highly variable genes, or use enrichment",
                "scores."), call. = FALSE)
        }
        m <- Matrix::t(.read_expr_matrix(ex, feats))
        i <- match(cells, rownames(m))
        if (anyNA(i)) {
            stop(fname, " ", nm, ": ", sum(is.na(i)), " cells have no ",
                "values", call. = FALSE)
        }
        m <- m[i, feats, drop = FALSE]
    }
    list(data = m, views = vw)
}

# dense W %*% M in blocks of columns
.misty_wprod <- function(W, M, block = 64L) {
    out <- matrix(0, nrow(W), ncol(M), dimnames = list(NULL, colnames(M)))
    for (s in seq(1L, max(1L, ncol(M)), by = block)) {
        k <- s:min(ncol(M), s + block - 1L)
        out[, k] <- as.matrix(W %*% M[, k, drop = FALSE])
    }
    out
}

.misty_view <- function(name, abbrev, data, own) {
    storage.mode(data) <- "double"
    list(name = name, abbrev = abbrev, data = data, own = own)
}

.misty_add_nbr <- function(views, W, M, tag, own, juxta = TRUE,
    para = TRUE) {
    if (juxta && !is.null(W$juxta)) {
        views[[length(views) + 1L]] <- .misty_view(
            paste0("juxtaview.", tag, W$juxta_radius),
            paste0("juxta.", tag, W$juxta_radius),
            .misty_wprod(W$juxta, M), own)
    }
    if (para) {
        for (l in names(W$para)) {
            views[[length(views) + 1L]] <- .misty_view(
                paste0("paraview.", tag, l), paste0("para.", tag, l),
                .misty_wprod(W$para[[l]], M), own)
        }
    }
    views
}

.misty_build_views <- function(gobject, spat_unit, feat_type, X, W,
    composition, extra_views, cells, fname, verbose) {
    views <- list(.misty_view("intraview", "intra", X, TRUE))
    views <- .misty_add_nbr(views, W, X, "", TRUE)
    if (!is.null(composition)) {
        lab <- .misty_cell_column(gobject, spat_unit, feat_type, composition,
            cells, fname)
        views <- .misty_add_nbr(views, W, .misty_onehot(lab), "composition.",
            FALSE)
    }
    for (nm in names(extra_views)) {
        src <- .misty_extra_source(gobject, spat_unit, feat_type,
            extra_views[[nm]], nm, cells, fname)
        if ("intra" %in% src$views) {
            views[[length(views) + 1L]] <- .misty_view(
                paste0("intraview.", nm), paste0("intra.", nm),
                as.matrix(src$data), FALSE)
        }
        views <- .misty_add_nbr(views, W, src$data, paste0(nm, "."), FALSE,
            juxta = "juxta" %in% src$views, para = "para" %in% src$views)
    }
    nms <- vapply(views, `[[`, "", "name")
    ab <- vapply(views, `[[`, "", "abbrev")
    if (anyDuplicated(nms) || anyDuplicated(ab)) {
        stop(fname, " duplicated view names: rename the extra views",
            call. = FALSE)
    }
    views
}

# mistyR draws its folds inside furrr (seed = TRUE), so under L'Ecuyer-CMRG
.misty_folds <- function(y, k, seed) {
    kind <- RNGkind()
    had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had) prev <- get(".Random.seed", envir = globalenv())
    on.exit({
        RNGkind(kind[1L], kind[2L], kind[3L])
        if (had) {
            assign(".Random.seed", prev, envir = globalenv())
        } else {
            rm(".Random.seed", envir = globalenv())
        }
    })
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    unname(caret::createFolds(y, k = k))
}

# lm() coefficients (NA where aliased) and optional t-test p-values
.misty_lm <- function(X, y, pvalues = FALSE) {
    X1 <- cbind(1, X)
    p <- ncol(X1)
    z <- .lm.fit(X1, y)
    r <- z$rank
    coef <- z$coefficients
    if (r < p) coef[(r + 1L):p] <- NA
    coef[z$pivot] <- coef
    out <- list(coef = coef)
    if (pvalues) {
        rdf <- length(y) - r
        pv <- rep(NA_real_, p)
        if (rdf > 0 && r > 0) {
            piv <- z$pivot[seq_len(r)]
            resvar <- sum(z$residuals^2) / rdf
            se <- sqrt(diag(chol2inv(z$qr[seq_len(r), seq_len(r),
                drop = FALSE])) * resvar)
            tv <- z$coefficients[seq_len(r)] / se
            pv[piv] <- 2 * stats::pt(abs(tv), rdf, lower.tail = FALSE)
        }
        out$pval <- pv[-1L]
    }
    out
}

.misty_lm_predict <- function(coef, X) {
    coef[is.na(coef)] <- 0
    drop(coef[1L] + X %*% coef[-1L])
}

# penalty of Cule and De Iorio (2012), test of Cule et al. (2011); no N x N U U'
.misty_ridge <- function(X, y) {
    n <- nrow(X)
    xm <- colMeans(X)
    ym <- mean(y)
    yc <- y - ym
    Xs <- sweep(X, 2L, xm)
    sc <- sqrt(colSums(Xs^2))
    s <- svd(sweep(Xs, 2L, sc, "/"))
    lam <- s$d^2
    uy <- drop(crossprod(s$u, yc))
    a <- uy / s$d
    share <- cumsum(lam) / sum(lam)
    r <- seq_len(max(1L, sum(share < 0.9) + any(share >= 0.9)))
    k_r <- r * (sum(yc^2) - cumsum(uy^2)[r]) / (n - r) / cumsum(a^2)[r]
    k_r <- k_r[seq_len(max(1L, sum(cumprod(is.finite(k_r)))))]
    df_r <- vapply(k_r, function(k) sum(lam^2 / (lam + k)^2), 1)
    k <- k_r[which.min(abs(df_r - seq_along(k_r)))]
    sh <- lam / (lam + k)
    b_s <- drop(s$v %*% (a * sh))
    beta <- b_s / sc
    rss <- sum((yc - s$u %*% (sh * uy))^2)
    s2 <- rss / (n - sum(lam * (lam + 2 * k) / (lam + k)^2))
    se <- sqrt(s2 * drop(s$v^2 %*% (lam / (lam + k)^2)))
    list(coef = c(ym - sum(beta * xm), beta),
        pval = 2 * (1 - stats::pnorm(abs(b_s / se))))
}

.misty_r2 <- function(pred, obs) {
    n <- sum(stats::complete.cases(pred))
    1 - (sum((obs - pred)^2) / ((n - 1) * stats::var(obs)))
}

.misty_rmse <- function(pred, obs) sqrt(mean((pred - obs)^2))

.misty_meta <- function(P, y, bypass, folds) {
    V <- ncol(P)
    use <- if (bypass) seq_len(V)[-1L] else seq_len(V)
    fit <- function(rows) {
        if (length(use) <= 1L) {
            .misty_lm(P[rows, use, drop = FALSE], y[rows],
                pvalues = TRUE)
        } else {
            .misty_ridge(P[rows, use, drop = FALSE], y[rows])
        }
    }
    all_rows <- seq_along(y)
    full <- fit(all_rows)
    b <- full$coef
    b[is.na(b)] <- 0
    p <- full$pval
    if (length(use) <= 1L) p[is.na(p)] <- 1
    coef <- if (bypass) c(b[1L], 0, b[-1L]) else b
    pval <- if (bypass) c(NA, 1, p) else c(NA, p)
    perf <- t(vapply(folds, function(f) {
        tr <- all_rows[-f]
        pi <- if (bypass) {
            rep(mean(y[tr]), length(f))
        } else {
            .misty_lm_predict(.misty_lm(P[tr, 1L, drop = FALSE],
                y[tr])$coef, P[f, 1L, drop = FALSE])
        }
        pm <- if (V == 1L) {
            pi
        } else {
            .misty_lm_predict(fit(tr)$coef, P[f, use, drop = FALSE])
        }
        c(intra.RMSE = .misty_rmse(pi, y[f]),
            intra.R2 = 100 * .misty_r2(pi, y[f]),
            multi.RMSE = .misty_rmse(pm, y[f]),
            multi.R2 = 100 * .misty_r2(pm, y[f]))
    }, numeric(4L)))
    negative <- any(perf < 0 | is.na(perf))
    perf[, c(2L, 4L)] <- pmax(perf[, c(2L, 4L)], 0, na.rm = TRUE)
    tt <- function(a, b, alt) {
        tryCatch(stats::t.test(a, b, alternative = alt)$p.value,
            error = function(e) 1)
    }
    summ <- c(colMeans(perf),
        p.RMSE = tt(perf[, 1L], perf[, 3L], "greater"),
        p.R2 = tt(perf[, 2L], perf[, 4L], "less"))
    summ[is.na(summ)] <- 1
    list(coef = coef, pval = pval, perf = summ, folds = perf,
        negative = negative)
}

.misty_view_model <- function(X, y, model, lm_folds, num_trees, mtry,
    sample_fraction, seed, threads) {
    if (model == "random_forest") {
        f <- ranger::ranger(x = X, y = y, num.trees = num_trees, mtry = mtry,
            sample.fraction = sample_fraction %null% 1,
            importance = "impurity", write.forest = FALSE, verbose = FALSE,
            num.threads = threads, seed = seed)
        return(list(pred = f$predictions, imp = f$variable.importance))
    }
    pred <- numeric(length(y))
    for (f in lm_folds) {
        pred[f] <- .misty_lm_predict(
            .misty_lm(X[-f, , drop = FALSE], y[-f])$coef,
            X[f, , drop = FALSE])
    }
    b <- .misty_lm(X, y)$coef[-1L]
    b[is.na(b)] <- 0
    names(b) <- colnames(X)
    list(pred = pred, imp = b)
}

.misty_fit_target <- function(target, Y, views, model, bypass, lm_folds,
    folds, num_trees, mtry, sample_fraction, seed, threads) {
    y <- Y[, target]
    n <- length(y)
    P <- matrix(0, n, length(views))
    imps <- vector("list", length(views))
    for (v in seq_along(views)) {
        X <- views[[v]]$data
        if (v == 1L && bypass) {
            X <- matrix(1, n, 1L, dimnames = list(NULL, ".novar"))
        } else if (views[[v]]$own) {
            X <- X[, colnames(X) != target, drop = FALSE]
        }
        if (!ncol(X)) {
            stop("[runMisty] view ", sQuote(views[[v]]$name),
                " has no predictors for target ", sQuote(target),
                call. = FALSE)
        }
        r <- .misty_view_model(X, y, model, lm_folds, num_trees, mtry,
            sample_fraction, seed, threads)
        P[, v] <- r$pred
        imps[[v]] <- r$imp
    }
    m <- .misty_meta(P, y, bypass, folds)
    m$imps <- imps
    m
}

.misty_fit_kth <- function(k, targets, meta_folds, ...) {
    .misty_fit_target(targets[k], folds = meta_folds[[k]], ...)
}

.misty_fit_sample <- function(X, views, targets, model, bypass, cv_folds,
    num_trees, mtry, sample_fraction, seed, threads, parallel = FALSE,
    verbose = FALSE, fname = "") {
    targets <- targets[apply(X[, targets, drop = FALSE], 2L, stats::sd) > 0]
    lm_folds <- if (model == "linear") {
        .misty_folds(seq_len(nrow(X)), cv_folds, seed)
    }
    meta_folds <- lapply(targets, function(t) {
        .misty_folds(X[, t], cv_folds, seed)
    })
    .misty_fit_targets(targets, X, views, model, bypass, lm_folds,
        meta_folds, num_trees, mtry, sample_fraction, seed, threads, parallel,
        verbose, fname)
}

.misty_fit_part <- function(part, ...) {
    .misty_fit_sample(part$X, part$views, ...)
}

.misty_fit_samples <- function(samples, X, views, targets, model, bypass,
    cv_folds, num_trees, mtry, sample_fraction, seed, cores, verbose,
    fname) {
    parts <- lapply(samples, function(r) {
        list(X = X[r, , drop = FALSE], views = lapply(views, function(v) {
            v$data <- v$data[r, , drop = FALSE]
            v
        }))
    })
    args <- list(targets = targets, model = model, bypass = bypass,
        cv_folds = cv_folds, num_trees = num_trees, mtry = mtry,
        sample_fraction = sample_fraction, seed = seed)
    par <- requireNamespace("future", quietly = TRUE) &&
        future::nbrOfWorkers() > 1
    if (par) {
        GiottoUtils::vmsg(.v = verbose, sprintf(
            "%s fitting %d samples on %d future workers", fname,
            length(parts), future::nbrOfWorkers()))
        old <- options(future.globals.maxSize = max(
            getOption("future.globals.maxSize", 0),
            as.numeric(utils::object.size(parts)) * 1.5))
        on.exit(options(old), add = TRUE)
        fits <- do.call(GiottoUtils::lapply_flex, c(list(X = parts,
            FUN = .misty_fit_part, method = "future", future.seed = TRUE,
            threads = 1L), args))
    } else {
        fits <- lapply(parts, function(p) {
            do.call(.misty_fit_part, c(list(p, threads = cores), args))
        })
    }
    names(fits) <- names(samples)
    fits
}

.misty_fit_notes <- function(fits, targets, verbose, fname) {
    n_neg <- sum(vapply(fits, function(f) {
        sum(vapply(f, `[[`, NA, "negative"))
    }, 1L))
    if (n_neg) {
        GiottoUtils::vmsg(.v = verbose, sprintf(paste(
            "%s %d target model(s) had negative or missing R2 in some folds,",
            "set to 0 as mistyR does"), fname, n_neg))
    }
    n_const <- sum(vapply(fits, function(f) {
        length(setdiff(targets, names(f)))
    }, 1L))
    if (n_const) {
        GiottoUtils::vmsg(.v = verbose, sprintf(paste(
            "%s %d target x sample pair(s) not fitted: the target is constant",
            "in the sample"), fname, n_const))
    }
    invisible(NULL)
}

.misty_fit_targets <- function(targets, X, views, model, bypass, lm_folds,
    meta_folds, num_trees, mtry, sample_fraction, seed, cores,
    parallel = FALSE, verbose = FALSE, fname = "") {
    if (!length(targets)) return(list())
    Y <- X[, targets, drop = FALSE]
    par <- parallel && requireNamespace("future", quietly = TRUE) &&
        future::nbrOfWorkers() > 1
    if (par) {
        GiottoUtils::vmsg(.v = verbose, sprintf(
            "%s fitting %d targets on %d future workers", fname,
            length(targets), future::nbrOfWorkers()))
        size <- (utils::object.size(views) + utils::object.size(Y)) * 1.5
        old <- options(future.globals.maxSize = max(
            getOption("future.globals.maxSize", 0), as.numeric(size)))
        on.exit(options(old), add = TRUE)
        fits <- GiottoUtils::lapply_flex(seq_along(targets),
            FUN = .misty_fit_kth, method = "future", future.seed = TRUE,
            targets = targets, meta_folds = meta_folds, Y = Y, views = views,
            model = model, bypass = bypass, lm_folds = lm_folds,
            num_trees = num_trees, mtry = mtry,
            sample_fraction = sample_fraction, seed = seed, threads = 1L)
        names(fits) <- targets
        return(fits)
    }
    fits <- vector("list", length(targets))
    for (k in seq_along(targets)) {
        t1 <- Sys.time()
        fits[[k]] <- .misty_fit_target(targets[k], Y, views, model, bypass,
            lm_folds, meta_folds[[k]], num_trees, mtry, sample_fraction, seed,
            cores)
        GiottoUtils::vmsg(.v = verbose, sprintf("%s %d/%d %s: %.1f s", fname,
            k, length(targets), targets[k],
            as.numeric(difftime(Sys.time(), t1, units = "secs"))))
    }
    names(fits) <- targets
    fits
}

.misty_write <- function(folder, targets, abbrev, fits) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    num <- function(x) paste(x, collapse = " ")
    writeLines(c(
        paste("target intercept", num(abbrev), "p.intercept",
            num(paste0("p.", abbrev))),
        vapply(targets, function(t) {
            paste(t, num(c(fits[[t]]$coef, fits[[t]]$pval)))
        }, "")
    ), file.path(folder, "coefficients.txt"))
    writeLines(c(
        "target intra.RMSE intra.R2 multi.RMSE multi.R2 p.RMSE p.R2",
        vapply(targets, function(t) paste(t, num(fits[[t]]$perf)), "")
    ), file.path(folder, "performance.txt"))
    for (t in targets) {
        for (v in seq_along(abbrev)) {
            imp <- fits[[t]]$imps[[v]]
            data.table::fwrite(
                data.table::data.table(target = names(imp),
                    imp = unname(imp)),
                file.path(folder, paste0("importances_", t, "_", abbrev[v],
                    ".txt")))
        }
    }
    invisible(folder)
}

# the list mistyR::collect_results() returns, over the samples
.misty_collect <- function(fits, abbrev) {
    measure <- value <- view <- Importance <- cv <- fraction <- NULL
    fits <- fits[lengths(fits) > 0L]
    tabs <- lapply(names(fits), function(s) {
        .misty_sample_tables(fits[[s]], abbrev, s)
    })
    bind <- function(x) data.table::rbindlist(lapply(tabs, `[[`, x))
    improvements <- bind("improvements")
    contributions <- bind("contributions")
    importances <- bind("importances")

    improvements.stats <- improvements[!startsWith(measure, "p."),
        list(mean = mean(value), sd = stats::sd(value)),
        keyby = c("target", "measure")]
    improvements.stats[, cv := sd / mean]
    data.table::setkey(improvements.stats, NULL)
    cs <- contributions[!startsWith(view, "p.") & view != "intercept",
        list(mean = mean(value)), keyby = c("target", "view")]
    cs[, fraction := abs(mean) / sum(abs(mean)), by = "target"]
    ps <- contributions[startsWith(view, "p.") & !grepl("intercept", view)]
    ps[, view := sub("^p\\.", "", view)]
    ps <- ps[, list(p.mean = mean(value), p.sd = stats::sd(value)),
        keyby = c("target", "view")]
    contributions.stats <- merge(cs, ps, by = c("target", "view"),
        sort = TRUE)
    data.table::setkey(contributions.stats, NULL)
    importances.aggregated <- importances[, list(
        Importance = mean(Importance), nsamples = .N),
        keyby = c("view", "Predictor", "Target")]
    data.table::setkey(importances.aggregated, NULL)

    out <- list(improvements = improvements, contributions = contributions,
        importances = importances, improvements.stats = improvements.stats,
        contributions.stats = contributions.stats,
        importances.aggregated = importances.aggregated)
    lapply(out, data.table::setDF)
}

.misty_sample_tables <- function(fits, abbrev, sample) {
    Importance <- intra.RMSE <- multi.RMSE <- intra.R2 <- multi.R2 <- NULL
    targets <- names(fits)
    perf <- data.table::rbindlist(lapply(targets, function(t) {
        data.table::as.data.table(as.list(fits[[t]]$perf))
    }))
    perf[, `:=`(target = targets, sample = sample)]
    perf[, `:=`(gain.RMSE = 100 * (intra.RMSE - multi.RMSE) / intra.RMSE,
        gain.R2 = multi.R2 - intra.R2)]
    m_cols <- c("intra.RMSE", "intra.R2", "multi.RMSE", "multi.R2", "p.RMSE",
        "p.R2", "gain.RMSE", "gain.R2")
    improvements <- data.table::data.table(
        target = rep(targets, each = length(m_cols)), sample = sample,
        measure = rep(m_cols, length(targets)),
        value = as.vector(t(as.matrix(perf[, m_cols, with = FALSE]))))

    v_cols <- c("intercept", abbrev, "p.intercept", paste0("p.", abbrev))
    contributions <- data.table::data.table(
        target = rep(targets, each = length(v_cols)), sample = sample,
        view = rep(v_cols, length(targets)),
        value = unlist(lapply(targets, function(t) {
            c(fits[[t]]$coef, fits[[t]]$pval)
        }), use.names = FALSE))

    tsort <- sort(targets, method = "radix")
    importances <- data.table::rbindlist(lapply(seq_along(abbrev), function(v) {
        imps <- lapply(tsort, function(t) fits[[t]]$imps[[v]])
        feats <- sort(unique(unlist(lapply(imps, names))), method = "radix")
        w <- vapply(seq_along(tsort), function(k) {
            imp <- imps[[k]][feats]
            p <- fits[[tsort[k]]]$pval[v + 1L]
            0 + scale(unname(imp))[, 1L] * (1 - p)
        }, numeric(length(feats)))
        w <- matrix(w, nrow = length(feats))
        data.table::data.table(sample = sample, view = abbrev[v],
            Predictor = rep(feats, each = length(tsort)),
            Target = rep(tsort, length(feats)),
            Importance = as.vector(t(w)))
    }))
    importances[is.nan(Importance), Importance := 0]
    list(improvements = improvements, contributions = contributions,
        importances = importances)
}

.misty_set_feat_meta <- function(gobject, spat_unit, feat_type, new_meta) {
    cols <- setdiff(colnames(new_meta), "feat_ID")
    fx <- GiottoClass::getFeatureMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, output = "featMetaObj", copy_obj = TRUE)
    old <- intersect(cols, colnames(fx[]))
    if (length(old)) {
        fx[] <- fx[][, setdiff(colnames(fx[]), old), with = FALSE]
        gobject <- GiottoClass::setFeatureMetadata(gobject, x = fx,
            verbose = FALSE)
    }
    GiottoClass::addFeatMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, new_metadata = new_meta, by_column = TRUE,
        column_feat_ID = "feat_ID")
}
