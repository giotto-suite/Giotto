# analyzeParam classes ####

#' @title Parameter Classes for Data Analysis Operations
#' @name analyze_param
#' @description
#' Parameter classes for use with [analyzeData()]. Each class encodes a
#' specific analysis method and its settings. Pass these to [analyzeData()]
#' directly or via the [analyzeParam()] factory function.
#'
#' These params produce computed values (scores, statistics). Any downstream
#' thresholding or selection is a separate step.
#'
#' **Stat params** (per-feature or per-cell summaries):
#' - `featStatsParam` — per-feature: nr_cells, perc_cells, mean_expr,
#'   mean_expr_det, total_expr
#' - `cellStatsParam` — per-cell: nr_feats, perc_feats, total_expr
#'
#' **COV-based score params** (coefficient of variation scores):
#' - `covGroupsParam` — COV z-score within expression-level bins
#' - `covLoessParam` — residual COV above a LOESS fit of COV ~ log(mean_expr)
#'
#' **Variance param**:
#' - `varParam` — per-feature variance of analytic Pearson residuals,
#'   computed from raw counts
#'
#' **Marker detection params** (one subclass of the virtual `markersParam`
#' per detection method):
#' - `scranMarkersParam` — pairwise group comparisons combined per group;
#'   see [markers_scran]. Built with [markersParam()].
#'
#' @param method character. One of `"feat_stats"`, `"cell_stats"`,
#'   `"cov_groups"`, `"cov_loess"`, `"var"`.
#' @param ... additional parameters passed to the specific param constructor.
#'   Use `$` on the returned object to inspect or modify individual params.
#' @returns an `analyzeParam`-inheriting object
#' @export
analyzeParam <- function(method, ...) {
    method <- match.arg(tolower(method),
        c("feat_stats", "cell_stats", "cov_groups", "cov_loess", "var")
    )
    switch(method,
        "feat_stats"  = .analyze_param_feat_stats(...),
        "cell_stats"  = .analyze_param_cell_stats(...),
        "cov_groups"  = .analyze_param_cov_groups(...),
        "cov_loess"   = .analyze_param_cov_loess(...),
        "var"         = .analyze_param_var(...)
    )
}

# VIRTUAL classes ####
#' @rdname analyze_param
#' @exportClass featStatsParam
setClass("featStatsParam", contains = "analyzeParam")
#' @rdname analyze_param
#' @exportClass cellStatsParam
setClass("cellStatsParam", contains = "analyzeParam")
#' @rdname analyze_param
#' @exportClass covGroupsParam
setClass("covGroupsParam", contains = "analyzeParam")
#' @rdname analyze_param
#' @exportClass covLoessParam
setClass("covLoessParam", contains = "analyzeParam")
#' @rdname analyze_param
#' @exportClass varParam
setClass("varParam", contains = "analyzeParam")

# constructors ####
.analyze_param_feat_stats <- function(...) {
    p <- new("featStatsParam", param = list(...))
    p$detection_threshold <- p$detection_threshold %null% 0
    p
}
.analyze_param_cell_stats <- function(...) {
    p <- new("cellStatsParam", param = list(...))
    p$detection_threshold <- p$detection_threshold %null% 0
    p
}
.analyze_param_cov_groups <- function(...) {
    p <- new("covGroupsParam", param = list(...))
    p$nr_expression_groups <- p$nr_expression_groups %null% 20
    p$detection_threshold <- p$detection_threshold %null% 0
    p$use_parallel <- p$use_parallel %null% FALSE
    p
}
.analyze_param_cov_loess <- function(...) {
    p <- new("covLoessParam", param = list(...))
    p$detection_threshold <- p$detection_threshold %null% 0
    p$use_parallel <- p$use_parallel %null% FALSE
    p
}
.analyze_param_var <- function(...) {
    p <- new("varParam", param = list(...))
    p$use_parallel <- p$use_parallel %null% FALSE
    # NB overdispersion for the Pearson residuals. 100 follows Lause/Kobak and
    # matches `normalizeGiotto(norm_methods = "pearson_resid")`, so both routes
    # to a residual evaluate one formula.
    p$theta <- p$theta %null% 100
    p
}



# Vectorized rowSds helper following flex function convention from GiottoClass
.rowSds_flex <- function(mymatrix, ...) {
    if (inherits(mymatrix, "IterableMatrix")) {
        v <- BPCells::matrix_stats(
            mymatrix, row_stats = "variance"
        )$row_stats["variance", ]
        return(sqrt(v))
    } else if (inherits(mymatrix, "DelayedArray")) {
        return(DelayedMatrixStats::rowSds(mymatrix, ...))
    } else if (inherits(mymatrix, "dgCMatrix")) {
        return(sparseMatrixStats::rowSds(mymatrix, ...))
    } else if (inherits(mymatrix, "Matrix")) {
        # For other Matrix types, use sparseMatrixStats
        return(sparseMatrixStats::rowSds(as(mymatrix, "dgCMatrix"), ...))
    } else if (inherits(mymatrix, "dbMatrix")) {
        # dbMatrix exports rowSds via MatrixGenerics
        return(dbMatrix::rowSds(mymatrix))
    } else {
        # Standard matrix - use matrixStats
        return(matrixStats::rowSds(as.matrix(mymatrix), ...))
    }
}


# Vectorized rowVars helper following flex function convention from GiottoClass
.rowVars_flex <- function(mymatrix, ...) {
    if (inherits(mymatrix, "DelayedArray")) {
        return(DelayedMatrixStats::rowVars(mymatrix, ...))
    } else if (inherits(mymatrix, "dgCMatrix")) {
        return(sparseMatrixStats::rowVars(mymatrix, ...))
    } else if (inherits(mymatrix, "Matrix")) {
        # For other Matrix types, use sparseMatrixStats
        return(sparseMatrixStats::rowVars(as(mymatrix, "dgCMatrix"), ...))
    } else if (inherits(mymatrix, "dbMatrix")) {
        # dbMatrix exports rowVars via MatrixGenerics
        return(dbMatrix::rowVars(mymatrix))
    } else {
        return(apply(mymatrix, 1, stats::var, ...))
    }
}


.calc_expr_general_stats <- function(expr_values, expression_threshold,
                                     calc_gini = TRUE) {
    # NSE vars
    gini <- NULL

    ## create data.table with relevant statistics ##
    # IterableMatrix fast path: a single BPCells::matrix_stats pass yields
    # mean, variance and nonzero — covering all four stats below when the
    # threshold is 0 (BPCells row_stats has no arbitrary-threshold count).
    if (inherits(expr_values, "IterableMatrix") &&
        expression_threshold == 0) {
        rs <- BPCells::matrix_stats(
            expr_values, row_stats = "variance"
        )$row_stats
        feat_in_cells_detected <- data.table::data.table(
            feats      = rownames(expr_values),
            nr_cells   = rs["nonzero", ],
            total_expr = rs["mean", ] * ncol(expr_values),
            mean_expr  = rs["mean", ],
            sd         = sqrt(rs["variance", ])
        )
    } else {
        feat_in_cells_detected <- data.table::data.table(
            feats = rownames(expr_values),
            nr_cells = rowSums_flex(expr_values > expression_threshold),
            total_expr = rowSums_flex(expr_values),
            mean_expr = rowMeans_flex(expr_values),
            sd = .rowSds_flex(expr_values)
        )
    }

    # calculate gini rowwise  (optional)
    if (isTRUE(calc_gini)) {
        gini_level <- unlist(apply(expr_values, MARGIN = 1, mygini_fun))
        feat_in_cells_detected[, gini := gini_level]
    }

    return(feat_in_cells_detected)
}


.calc_expr_cov_stats <- function(expr_values, expression_threshold,
                                  calc_gini = TRUE) {
    # NSE vars
    cov <- sd <- mean_expr <- NULL

    # get general expression statistics and gini data.table
    feat_in_cells_detected <- .calc_expr_general_stats(
        expr_values, expression_threshold, calc_gini = calc_gini
    )

    # calculate cov using sd and mean_expr from general stats DT
    feat_in_cells_detected[, cov := (sd / mean_expr)]

    return(feat_in_cells_detected)
}


.calc_expr_cov_stats_parallel <- function(expr_values,
    expression_threshold,
    calc_gini = TRUE,
    cores = GiottoUtils::determine_cores()) {
    # NSE vars
    cov <- sd <- mean_expr <- NULL

    # setup chunk rows to use for each parallel based on number of cores
    chunk_rows <- seq(nrow(expr_values)) %>%
        split(., cut(., cores))

    # params to pass into the future_lapply
    fparams <- list(
        calc_fun = .calc_expr_general_stats,
        expression_threshold = expression_threshold,
        calc_gini = calc_gini
    )

    # parallelized calculation of general stats
    chunk_stats_dt_list <- lapply_flex(
        chunk_rows,
        function(r_idx, fparams) {
            fparams$calc_fun(
                expr_values = expr_values[r_idx, ],
                expression_threshold = fparams$expression_threshold,
                calc_gini = fparams$calc_gini
            )
        },
        fparams = fparams,
        cores = cores,
        future.seed = TRUE
    )

    # combine stats tables
    feat_in_cells_detected <- data.table::rbindlist(chunk_stats_dt_list)

    # calculate cov using sd and mean_expr from combined general stats DT
    feat_in_cells_detected[, cov := (sd / mean_expr)]

    return(feat_in_cells_detected)
}





#' @title calculateHVF
#' @name calculateHVF
#' @description compute highly variable features
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param method method to calculate highly variable features
#' @param reverse_log_scale reverse log-scale of expression values
#' (default = FALSE)
#' @param logbase if `reverse_log_scale` is TRUE, which log base was used?
#' @param expression_threshold expression threshold to consider a gene detected
#' @param nr_expression_groups (cov_groups) number of expression groups for
#' cov_groups
#' @param zscore_threshold (cov_groups) zscore to select hvg for cov_groups
#' @param HVFname name for highly variable features in cell metadata
#' @param difference_in_cov (cov_loess) minimum difference in coefficient of
#' variance required
#' @param var_threshold (var_p_resid) noise floor on the residual variance.
#' Under analytic Pearson residuals a feature with no biological signal has
#' variance close to 1 by construction, so the default `1` drops features that
#' are *less* variable than Poisson noise. It is a guard, not a selector: on
#' real data it keeps most features (94.5% of a Stereo-seq transcriptome, 82%
#' of a Xenium panel). Raise it to select on variance alone. `NULL` disables it.
#' @param n_top_feats maximum number of features to keep, taken from the top of
#' that method's own ranking. Default `2000`. Applies to every method.
#' `NULL` disables it.
#' @param var_number `r lifecycle::badge("deprecated")` use `n_top_feats`.
#' @details
#' # Feature selection
#' Every method applies two constraints and the more restrictive one decides:
#' its own threshold (`zscore_threshold` for `cov_groups`,
#' `difference_in_cov` for `cov_loess`, `var_threshold` for `var_p_resid`)
#' and the shared `n_top_feats` count. Setting either to `NULL` falls back to
#' the other.
#'
#' The two cover different failure modes. A count keeps the selection a
#' predictable size across samples, but pads the list when a dataset has few
#' genuinely variable features, and cannot bind at all on a panel smaller than
#' `n_top_feats`. A threshold is scale-free but its yield swings with the
#' data.
#'
#' Each method ranks on its own score — residual variance, within-bin COV
#' z-score, or COV above the LOESS fit. For `cov_groups` the z-score is
#' standardised within expression bins, so a global ranking stays balanced
#' across expression levels: on a Stereo-seq transcriptome a top-2000 cut drew
#' 9.6-12.4% from each of 19 of 20 bins, the exception being the
#' highest-expression bin, where COV is compressed by construction.
#' Ties resolve by `order()`, so the selected count is exact.
#' @param theta (var_p_resid) negative-binomial overdispersion used when
#' forming the Pearson residuals. Default 100, following Lause/Kobak and
#' matching `normalizeGiotto(norm_methods = "pearson_resid")`. Larger values
#' approach a Poisson model.
#' @param random_subset random subset to perform HVF detection on.
#' Passing `NULL` runs HVF on all cells.
#' @param set_seed logical. whether to set a seed when random_subset is used
#' @param seed_number seed number to use when random_subset is used
#' @param show_plot show plot
#' @param return_plot return ggplot object (overridden by `return_gobject`)
#' @param save_plot logical. directly save the plot
#' @param save_param list of saving parameters from
#' [GiottoVisuals::all_plots_save_function()]
#' @param default_save_name default save name for saving, don't change, change
#' save_name in save_param
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param calc_gini logical. Whether to calculate Gini index for each feature.
#' Set to FALSE for performance with large datasets or dbMatrix objects.
#' @param verbose be verbose
#' @returns giotto object highly variable features appended to feature metadata
#' (`fDataDT()`)
#' @details
#' Currently we provide 2 ways to calculate highly variable genes:
#'
#' \strong{1. high coeff of variance (COV) within groups: } \cr
#' First genes are binned (\emph{nr_expression_groups}) into average expression
#' groups and the COV for each feature is converted into a z-score within each
#' bin. Features with a z-score higher than the threshold
#' (\emph{zscore_threshold}) are considered highly variable.  \cr
#'
#' \strong{2. high COV based on loess regression prediction: } \cr
#' A predicted COV is calculated for each feature using loess regression
#' (COV~log(mean expression))
#' Features that show a higher than predicted COV (\emph{difference_in_cov})
#' are considered highly variable. \cr
#'
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' calculateHVF(g)
#' @export
calculateHVF <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        method = c("cov_groups", "cov_loess", "var_p_resid"),
        reverse_log_scale = FALSE,
        logbase = 2,
        expression_threshold = 0,
        nr_expression_groups = 20,
        zscore_threshold = 1.5,
        HVFname = "hvf",
        difference_in_cov = 0.1,
        var_threshold = 1,
        n_top_feats = 2000,
        theta = 100,
        var_number = deprecated(),
        random_subset = NULL,
        set_seed = TRUE,
        seed_number = 1234,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "HVFplot",
        return_gobject = TRUE,
        calc_gini = FALSE,
        verbose = TRUE) {
    # NSE vars
    selected <- feats <- var <- NULL

    # `var_number` was var_p_resid-only; the top-N cut now applies to every
    # method, so it is exposed under a name that does not imply variance.
    n_top_feats <- GiottoUtils::deprecate_param(
        var_number, n_top_feats, fun = "calculateHVF", when = "4.2.4"
    )

    # determine whether to use parallel functions
    # Do not use future if future packages are not installed
    # Do not use future if plan is "sequential"
    has_future <- requireNamespace("future.apply", quietly = TRUE) &&
        requireNamespace("future", quietly = TRUE)
    use_parallel <- ifelse(has_future,
        !("sequential" %in% class(future::plan())),
        FALSE
    )

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # Resolved here rather than further down, because which expression values
    # get fetched depends on it.
    method <- match.arg(
        method,
        choices = c("cov_groups", "cov_loess", "var_p_resid")
    )

    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )

    # `var_p_resid` builds Pearson residuals itself, and they are only defined
    # on counts, so it always reads raw regardless of `expression_values`.
    # This used to require running
    # `normalizeGiotto(norm_methods = "pearson_resid")` first and then asking
    # for the slot it wrote to -- "scaled", not the "normalized" default --
    # a three-way agreement nothing checked, and which silently returned the
    # variance of library-normalized values when it was not met.
    if (method == "var_p_resid") {
        if (!missing(expression_values)) {
            warning(wrap_txt(
                "[calculateHVF] expression_values =", paste0("'", values, "'"),
                "is ignored for method = 'var_p_resid': Pearson residuals are",
                "computed from raw counts. Use normalizeGiotto(norm_methods =",
                "'pearson_resid') if you want residuals as a stored matrix."
            ), call. = FALSE)
        }
        values <- "raw"
    }

    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )

    # not advised
    if (isTRUE(reverse_log_scale)) {
        expr_values <- (logbase^expr_values) - 1
    }

    # create a random subset if random_subset is not NULL
    if (!is.null(random_subset)) {
        if (isTRUE(set_seed)) GiottoUtils::local_seed(seed = seed_number)

        random_selection <- sort(sample(
            seq_len(ncol(expr_values)), random_subset
        ))
        expr_values <- expr_values[, random_selection]

        if (isTRUE(set_seed)) GiottoUtils::random_seed()
    }



    # print, return and save parameters
    show_plot <- ifelse(is.na(show_plot),
        instructions(gobject, "show_plot"),
        show_plot
    )
    save_plot <- ifelse(is.na(save_plot),
        instructions(gobject, "save_plot"),
        save_plot
    )
    return_plot <- ifelse(is.na(return_plot),
        instructions(gobject, "return_plot"),
        return_plot
    )


    # `method` is already resolved above, where it selects the input values.

    # Stats compute is dispatched on the expression backend via
    # analyzeData(x, analyzeParam(method, ...)). Streaming backends
    # (parquetExprStore in GiottoDisk) provide their own setMethod for the
    # same generic and inherit this user-facing API.
    want_plot <- any(isTRUE(show_plot), isTRUE(return_plot), isTRUE(save_plot))
    cov_diff <- cov_group_zscore <- expr_groups <- mean_expr <- pred_cov_feats <-
        rank <- selected <- var <- cov <- NULL  # NSE bindings

    feat_in_cells_detected <- switch(method,
        "var_p_resid" = {
            dt <- analyzeData(
                expr_values,
                analyzeParam("var", use_parallel = use_parallel,
                             theta = theta)
            )
            keep <- if (is.null(var_threshold)) {
                rep(TRUE, nrow(dt))
            } else {
                dt$var >= var_threshold
            }
            dt[, selected := .hvf_select(keep, dt$var, n_top_feats)]
            dt
        },
        "cov_groups" = {
            dt <- analyzeData(
                expr_values,
                analyzeParam("cov_groups",
                    nr_expression_groups = nr_expression_groups,
                    detection_threshold = expression_threshold,
                    use_parallel = use_parallel)
            )
            dt[, selected := .hvf_select(
                cov_group_zscore > zscore_threshold,
                cov_group_zscore, n_top_feats
            )]
            dt
        },
        "cov_loess" = {
            dt <- analyzeData(
                expr_values,
                analyzeParam("cov_loess",
                    detection_threshold = expression_threshold,
                    use_parallel = use_parallel)
            )
            dt[, selected := .hvf_select(
                cov_diff > difference_in_cov, cov_diff, n_top_feats
            )]
            dt
        }
    )

    # Plot generation. The analyzeData methods drop intermediate columns
    # (expr_groups, pred_cov_feats, rank) used only for plotting; re-add
    # them locally so the existing plot helpers can be used unchanged.
    pl <- NULL
    if (want_plot) {
        if (method == "var_p_resid") {
            # ggplot holds a reference to the data, and data.table's `:=`
            # edits in place -- dropping `rank` after the fact would strip it
            # back out of the plot. Hand the plot its own copy.
            feat_in_cells_detected[, rank := seq_len(.N)]
            pl <- .create_calc_var_hvf_plot(
                data.table::copy(feat_in_cells_detected),
                var_threshold = var_threshold,
                n_top_feats = n_top_feats
            )
            feat_in_cells_detected[, rank := NULL]
        } else if (method == "cov_groups") {
            prob_seq <- seq(0, 1, 1 / nr_expression_groups)
            prob_seq[length(prob_seq)] <- 1
            breaks <- stats::quantile(feat_in_cells_detected$mean_expr,
                                       probs = prob_seq)
            if (any(duplicated(breaks))) {
                v <- feat_in_cells_detected$mean_expr
                breaks <- stats::quantile(v[v > 0], probs = prob_seq)
                breaks[[1]] <- 0
            }
            feat_in_cells_detected[, expr_groups := cut(
                mean_expr, breaks = breaks,
                labels = paste0("group_", seq_len(nr_expression_groups)),
                include.lowest = TRUE
            )]
            pl <- .create_cov_group_hvf_plot(
                feat_in_cells_detected, nr_expression_groups
            )
            feat_in_cells_detected[, expr_groups := NULL]
        } else if (method == "cov_loess") {
            loess_model <- stats::loess(cov ~ log(mean_expr),
                                         data = feat_in_cells_detected)
            feat_in_cells_detected[, pred_cov_feats := stats::predict(
                loess_model, newdata = feat_in_cells_detected
            )]
            pl <- .create_cov_loess_hvf_plot(
                feat_in_cells_detected, difference_in_cov, var_col = "cov"
            )
            feat_in_cells_detected[, pred_cov_feats := NULL]
        }
    }




    ## print plot
    if (isTRUE(show_plot)) {
        print(pl)
    }

    ## save plot
    if (isTRUE(save_plot)) {
        do.call(
            GiottoVisuals::all_plots_save_function,
            c(list(
                gobject = gobject, plot_object = pl,
                default_save_name = default_save_name
            ), save_param)
        )
    }

    ## return plot
    if (isTRUE(return_plot)) {
        if (isTRUE(return_gobject)) {
            message("return_plot = TRUE and return_gobject = TRUE \n
                    plot will not be returned to object, but can still be
                    saved with save_plot = TRUE or manually")
        } else {
            return(pl)
        }
    }


    if (isTRUE(return_gobject)) {
        # add HVG metadata to feat_metadata
        feat_metadata <- getFeatureMetadata(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "featMetaObj",
            copy_obj = TRUE
        )

        column_names_feat_metadata <- colnames(feat_metadata[])

        if (HVFname %in% column_names_feat_metadata) {
            vmsg(
                .v = verbose, HVFname,
                " has already been used, will be overwritten"
            )
            feat_metadata[][, eval(HVFname) := NULL]

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- setFeatureMetadata(gobject,
                x = feat_metadata,
                verbose = FALSE,
                initialize = FALSE
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }

        if (method == "var_p_resid") {
            HVGfeats <- feat_in_cells_detected[, .(feats, var, selected)]
            data.table::setnames(HVGfeats, "selected", HVFname)
        } else {
            HVGfeats <- feat_in_cells_detected[, .(feats, selected)]
            data.table::setnames(HVGfeats, "selected", HVFname)
        }


        gobject <- addFeatMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = HVGfeats,
            by_column = TRUE,
            column_feat_ID = "feats"
        )

        ## update parameters used ##
        gobject <- update_giotto_params(gobject, description = "_hvf")

        return(gobject)
    } else {
        return(feat_in_cells_detected)
    }
}








# plot generation ####
.create_cov_group_hvf_plot <- function(feat_in_cells_detected,
    nr_expression_groups) {
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic() +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = 14),
            axis.text = ggplot2::element_text(size = 12)
        )
    pl <- pl + ggplot2::geom_point(
        data = feat_in_cells_detected,
        GiottoVisuals::aes_string2(
            x = "mean_expr", y = "cov", color = "selected")
    )
    pl <- pl + ggplot2::scale_color_manual(
        values = c(no = "lightgrey", yes = "orange"),
        guide = ggplot2::guide_legend(
            title = "HVF",
            override.aes = list(size = 5)
        )
    )
    pl <- pl + ggplot2::facet_wrap(
        ~expr_groups,
        ncol = nr_expression_groups, scales = "free_x"
    )
    pl <- pl + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 4)
    )
    pl <- pl + ggplot2::labs(x = "expression groups", y = "cov")
    pl
}


.create_cov_loess_hvf_plot <- function(feat_in_cells_detected,
    difference_in_cov, var_col) {
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic() +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = 14),
            axis.text = ggplot2::element_text(size = 12)
        )
    pl <- pl + ggplot2::geom_point(
        data = feat_in_cells_detected,
        GiottoVisuals::aes_string2(
            x = "log(mean_expr)", y = var_col, color = "selected")
    )
    pl <- pl + ggplot2::geom_line(
        data = feat_in_cells_detected,
        GiottoVisuals::aes_string2(x = "log(mean_expr)", y = "pred_cov_feats"),
        color = "blue"
    )
    hvg_line <- paste0("pred_cov_feats+", difference_in_cov)
    pl <- pl + ggplot2::geom_line(
        data = feat_in_cells_detected,
        GiottoVisuals::aes_string2(x = "log(mean_expr)", y = hvg_line), 
        linetype = 2
    )
    pl <- pl + ggplot2::labs(x = "log(mean expression)", y = var_col)
    pl <- pl + ggplot2::scale_color_manual(
        values = c(no = "lightgrey", yes = "orange"),
        guide = ggplot2::guide_legend(
            title = "HVF",
            override.aes = list(size = 5)
        )
    )
    pl
}


# Diagnostic for method = "var_p_resid".
#
# Two things make this readable that a bare rank-vs-variance scatter does not:
# a reference line at var = 1 -- the value a feature with no biological signal
# takes under analytic Pearson residuals, so it is the yardstick the threshold
# is relative to -- and a log y axis, because residual variances are
# heavy-tailed and the elbow is invisible on a linear scale.
#
# When `mean_expr` is available a second panel plots it against the residual
# variance. Decoupling variance from expression level is the whole point of
# Pearson residuals, so that panel is what shows whether the method did its
# job: selected features should not bunch at high mean.
# Apply the shared top-N cut on top of a method's own threshold, and return
# the "yes"/"no" column.
#
# `keep` is the method's threshold verdict; `score` is the quantity it ranks
# on (residual variance, within-bin COV z-score, COV above the LOESS fit).
# Both constraints apply and the more restrictive one wins. They cover
# different failure modes: a count keeps the selection a predictable size
# across samples but pads the list when a dataset has few variable features --
# and cannot bind at all on a panel smaller than `n_top`; a threshold is
# scale-free but its yield swings with the data (the old var_p_resid default
# of 1.5 took 49.6% of a Stereo-seq transcriptome and 39.2% of a Xenium panel).
#
# The rank is computed here rather than taken from row order, because the
# `analyzeData()` methods do not agree on it: cov_loess and var return sorted,
# cov_groups does not. Ties resolve by `order()`, so the count is exact.
.hvf_select <- function(keep, score, n_top) {
    if (!is.null(n_top) && is.numeric(n_top)) {
        rank <- integer(length(score))
        rank[order(-score)] <- seq_along(score)
        keep <- keep & (rank <= min(as.integer(n_top), length(score)))
    }
    ifelse(keep, "yes", "no")
}


.create_calc_var_hvf_plot <- function(dt_res, var_threshold = NULL,
                                       n_top_feats = NULL) {
    n_sel <- sum(dt_res$selected == "yes")
    n_all <- nrow(dt_res)

    # Both constraints apply; say which one actually decided, since that is
    # the thing a user needs to know before touching either.
    bound <- if (!is.null(n_top_feats) && n_sel == min(n_top_feats, n_all)) {
        sprintf("count-limited (top %s)", format(min(n_top_feats, n_all),
                                                  big.mark = ","))
    } else if (!is.null(var_threshold)) {
        sprintf("threshold-limited (var >= %g)", var_threshold)
    } else {
        "unconstrained"
    }
    sub <- sprintf("%s of %s features selected -- %s",
                   format(n_sel, big.mark = ","),
                   format(n_all, big.mark = ","), bound)

    .base <- function(pl) {
        pl <- pl + ggplot2::geom_hline(
            yintercept = 1, linetype = "solid", colour = "steelblue"
        )
        if (!is.null(var_threshold)) {
            pl <- pl + ggplot2::geom_hline(
                yintercept = var_threshold, linetype = "dashed",
                colour = "black"
            )
        }
        pl + ggplot2::scale_y_log10() +
            ggplot2::scale_color_manual(
                values = c(no = "lightgrey", yes = "orange"),
                guide = ggplot2::guide_legend(
                    title = "HVF", override.aes = list(size = 5)
                )
            ) +
            ggplot2::theme_classic() +
            ggplot2::theme(
                axis.title = ggplot2::element_text(size = 14),
                axis.text = ggplot2::element_text(size = 12)
            )
    }

    pl_rank <- .base(
        ggplot2::ggplot() + ggplot2::geom_point(
            data = dt_res,
            GiottoVisuals::aes_string2(
                x = "rank", y = "var", color = "selected"
            )
        ) + ggplot2::scale_x_reverse()
    )
    # The count cut is a rank, so it only has a place on this panel.
    if (!is.null(n_top_feats) && n_top_feats < n_all) {
        pl_rank <- pl_rank + ggplot2::geom_vline(
            xintercept = n_top_feats, linetype = "dotted", colour = "black"
        )
    }
    pl_rank <- pl_rank + ggplot2::labs(
        x = "feature rank", y = "Pearson residual variance",
        title = "HVF by residual variance", subtitle = sub,
        caption = paste("solid: var = 1, the no-signal expectation;",
                        "dashed: var_threshold; dotted: n_top_feats")
    )

    if (!"mean_expr" %in% names(dt_res)) return(pl_rank)

    pl_mean <- .base(
        ggplot2::ggplot() + ggplot2::geom_point(
            data = dt_res,
            GiottoVisuals::aes_string2(
                x = "mean_expr", y = "var", color = "selected"
            )
        ) + ggplot2::scale_x_log10()
    ) + ggplot2::labs(
        x = "mean expression", y = "Pearson residual variance",
        title = "Residual variance vs expression level", subtitle = sub,
        caption = "selection should not track mean expression"
    )

    if (requireNamespace("patchwork", quietly = TRUE)) {
        return(pl_rank + pl_mean)
    }
    # patchwork is optional; the rank view is the primary one
    pl_rank
}


# analyzeData methods ####

#' @name analyzeData
#' @title Data Analysis via Parameter Dispatch
#' @description
#' Compute statistics or scores from matrix-type data. `analyzeData()` is a
#' generic that dispatches on both `x` (the data) and `param` (the analysis
#' operation). Methods return a `data.table` of computed values; any downstream
#' thresholding or selection is a separate step.
#' @param x data to analyze
#' @param param an [analyzeParam-class] inheriting object defining the analysis
#' operation and its settings
#' @param \dots additional params passed to specific methods
#' @returns a `data.table` of computed values
NULL

# exprObj base dispatch
#' @rdname analyzeData
setMethod("analyzeData",
    signature(x = "exprObj", param = "analyzeParam"),
    function(x, param, ...) {
        analyzeData(x[], param, ...)
    }
)

# * featStatsParam ####
#' @rdname analyzeData
#' @param groups optional vector of group assignments, one per column of `x`,
#' `NA` to exclude. When supplied, the statistics are taken per
#' (feature, group) instead of over every cell, and the result gains `group`
#' and `n_cells` columns.
#' @param stats optional character vector of accumulators to compute, any of
#' `"sum"`, `"sumsq"`, `"nnz"`, `"sum_det"`. Grouped path only; emitted columns
#' are whichever the requested accumulators support.
setMethod("analyzeData",
    signature(x = "allMatrix", param = "featStatsParam"),
    function(x, param, ..., groups = NULL, stats = NULL) {
        det_thresh <- param$detection_threshold

        if (!is.null(groups)) {
            return(.feat_stats_grouped(
                x, det_thresh, groups,
                stats = stats %null% c("sum", "sumsq", "nnz", "sum_det")
            ))
        }
        if (!is.null(stats)) {
            stop("[feat_stats] `stats` selection is only available on the ",
                "grouped path (pass `groups`). The ungrouped verb has a ",
                "fixed column contract.", call. = FALSE)
        }

        mean_expr_det <- NULL
        n_detected <- rowSums_flex(x > det_thresh)
        feat_stats <- data.table::data.table(
            feats      = rownames(x),
            nr_cells   = n_detected,
            perc_cells = (n_detected / ncol(x)) * 100,
            total_expr = rowSums_flex(x),
            mean_expr  = rowMeans_flex(x)
        )
        feat_stats[, mean_expr_det := .mean_expr_det_test(
            x, detection_threshold = det_thresh
        )]
        feat_stats
    }
)

# Per-(feature, group) statistics: the same accumulators, partitioned by a
# per-cell grouping instead of taken over every cell.
#
# The contract matches GiottoDisk's streaming implementation exactly, because
# the whole point is that a caller cannot tell the backends apart -- same
# columns, same zero-filled feats x groups cross product, same `adr/0009`
# threshold semantics where the detection threshold gates `nr_cells` but never
# the sums.
#
# Emitting the complete cross product matters: a feature with no expression in
# a group has mean 0 over that group's cells, not a missing row. Gini is taken
# over the length-G vector per feature, so a dropped row would silently change
# the coefficient.
.feat_stats_grouped <- function(x, thr, groups,
    stats = c("sum", "sumsq", "nnz", "sum_det")) {
    stats <- match.arg(stats, several.ok = TRUE)

    n_cells <- ncol(x)
    if (length(groups) != n_cells) {
        stop("[feat_stats] `groups` must have one entry per cell (",
            n_cells, "), got ", length(groups), ".", call. = FALSE)
    }

    # `droplevels` so an unused level cannot surface as a group of zero cells
    g <- droplevels(if (is.factor(groups)) groups else factor(groups))
    lvls <- levels(g)
    if (length(lvls) < 1L) {
        stop("[feat_stats] `groups` has no non-empty levels.", call. = FALSE)
    }

    n_feats <- nrow(x)
    nk <- as.numeric(tabulate(as.integer(g), nbins = length(lvls)))
    names(nk) <- lvls

    acc <- lapply(stats, function(nm) numeric(0))
    names(acc) <- stats
    per <- lapply(lvls, function(k) {
        sub <- x[, which(g == k), drop = FALSE]
        det <- if (any(c("nnz", "sum_det") %in% stats)) sub > thr
        list(
            # `rowMeans_flex` rather than sum/n: this is what
            # `create_average_DT()` used, so the rewired gini path stays
            # bit-identical rather than merely equal to tolerance.
            mean = if ("sum" %in% stats) rowMeans_flex(sub),
            sum = if ("sum" %in% stats) rowSums_flex(sub),
            sumsq = if ("sumsq" %in% stats) rowSums_flex(sub * sub),
            nnz = if ("nnz" %in% stats) rowSums_flex(det),
            sum_det = if ("sum_det" %in% stats) rowSums_flex(sub * det)
        )
    })

    # groups slowest, feats cycling within -- the order `matrix` unrolls in,
    # matching what the streaming backend emits
    nn <- rep(nk, each = n_feats)
    out <- data.table::data.table(
        feats = rep(rownames(x), times = length(lvls)),
        group = rep(lvls, each = n_feats),
        n_cells = nn
    )
    pull <- function(nm) as.numeric(unlist(lapply(per, `[[`, nm)))

    if ("sum" %in% stats) {
        gene_sum <- pull("sum")
        out[, "total_expr" := gene_sum]
        out[, "mean_expr" := pull("mean")]

        if ("sumsq" %in% stats) {
            gene_sumsq <- pull("sumsq")
            out[, "sumsq" := gene_sumsq]
            # clamped: the subtraction can go slightly negative when the mean
            # dominates the spread, and a negative variance would surface as
            # an NaN standard deviation
            gene_var <- ifelse(nn > 1,
                pmax((gene_sumsq - gene_sum * gene_sum / nn) / (nn - 1), 0), 0
            )
            out[, "sd" := sqrt(gene_var)]
        }
    }
    if ("nnz" %in% stats) {
        gene_nnz <- pull("nnz")
        out[, "nr_cells" := as.integer(gene_nnz)]
        out[, "perc_cells" := ifelse(nn > 0, gene_nnz / nn * 100, NaN)]

        if ("sum_det" %in% stats) {
            out[, "mean_expr_det" := ifelse(
                gene_nnz > 0, pull("sum_det") / gene_nnz, NaN
            )]
        }
    }

    out[]
}

# * cellStatsParam ####
#' @rdname analyzeData
setMethod("analyzeData",
    signature(x = "allMatrix", param = "cellStatsParam"),
    function(x, param, ...) {
        det_thresh <- param$detection_threshold
        n_detected <- colSums_flex(x > det_thresh)
        data.table::data.table(
            cells      = colnames(x),
            nr_feats   = n_detected,
            perc_feats = (n_detected / nrow(x)) * 100,
            total_expr = colSums_flex(x)
        )
    }
)

# * covGroupsParam ####
#' @rdname analyzeData
setMethod("analyzeData",
    signature(x = "allMatrix", param = "covGroupsParam"),
    function(x, param, ...) {
        cov_group_zscore <- cov <- expr_groups <- NULL
        nr_groups <- param$nr_expression_groups
        det_thresh <- param$detection_threshold

        # BPCells already streams in C++ with internal threading; R-level
        # chunking via future_lapply re-opens the file per worker and is a
        # net loss for IterableMatrix.
        calc_fun <- if (isTRUE(param$use_parallel) &&
                        !inherits(x, "IterableMatrix")) {
            .calc_expr_cov_stats_parallel
        } else {
            .calc_expr_cov_stats
        }
        dt <- calc_fun(x, expression_threshold = det_thresh, calc_gini = FALSE)
        dt <- dt[nr_cells > 0]

        prob_sequence <- seq(0, 1, 1 / nr_groups)
        prob_sequence[length(prob_sequence)] <- 1
        breaks <- stats::quantile(dt$mean_expr, probs = prob_sequence)
        if (any(duplicated(breaks))) {
            v <- dt$mean_expr
            breaks <- stats::quantile(v[v > 0], probs = prob_sequence)
            breaks[[1]] <- 0
        }
        dt[, expr_groups := cut(
            mean_expr, breaks = breaks,
            labels = paste0("group_", seq_len(nr_groups)),
            include.lowest = TRUE
        )]
        dt[, cov_group_zscore := scale(cov), by = expr_groups]
        dt[, expr_groups := NULL]
        dt
    }
)

# * covLoessParam ####
#' @rdname analyzeData
setMethod("analyzeData",
    signature(x = "allMatrix", param = "covLoessParam"),
    function(x, param, ...) {
        pred_cov <- cov_diff <- NULL
        det_thresh <- param$detection_threshold

        # BPCells already streams in C++ with internal threading; R-level
        # chunking via future_lapply re-opens the file per worker and is a
        # net loss for IterableMatrix.
        calc_fun <- if (isTRUE(param$use_parallel) &&
                        !inherits(x, "IterableMatrix")) {
            .calc_expr_cov_stats_parallel
        } else {
            .calc_expr_cov_stats
        }
        dt <- calc_fun(x, expression_threshold = det_thresh, calc_gini = FALSE)
        dt <- dt[nr_cells > 0]

        loess_fit <- stats::loess(cov ~ log(mean_expr), data = dt)
        dt[, pred_cov  := stats::predict(loess_fit, newdata = dt)]
        dt[, cov_diff  := cov - pred_cov]
        dt[, pred_cov  := NULL]
        data.table::setorder(dt, -cov_diff)
        dt
    }
)

# * varParam ####
#' @rdname analyzeData
setMethod("analyzeData",
    signature(x = "allMatrix", param = "varParam"),
    function(x, param, ...) {
        # Analytic Pearson residual variance (Lause/Kobak/Berens), computed
        # here from raw counts rather than expecting the caller to have run
        # `normalizeGiotto(norm_methods = "pearson_resid")` first. Previously
        # this was a plain `var` of whatever matrix it was handed, which was
        # the residual criterion only if that prior step had happened and had
        # been read back from the right slot -- it silently was not.
        #
        # `.prnorm()` is the same helper the normalization uses, so the two
        # routes cannot drift apart.
        theta <- param$theta %null% 100

        if (inherits(x, "IterableMatrix")) {
            # BPCells has no residual path; rowVars on the stored values is
            # not the same statistic, so say so rather than return it.
            stop("[analyzeData(varParam)] Pearson residual variance is not ",
                 "implemented for IterableMatrix. Use method = 'cov_loess' ",
                 "or 'cov_groups', or a parquetExprStore backend.",
                 call. = FALSE)
        }

        z <- .prnorm(x = x, theta = theta)
        scores <- .rowVars_flex(z)
        names(scores) <- rownames(x)
        mu <- rowMeans_flex(x)

        dt <- data.table::data.table(
            feats = names(scores),
            var = as.numeric(scores),
            mean_expr = as.numeric(mu[names(scores)])
        )
        data.table::setorder(dt, -var)
        dt
    }
)

