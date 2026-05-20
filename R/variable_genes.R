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
#' - `varParam` — per-feature variance on a scaled matrix
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
#' @param var_threshold (var_p_resid) variance threshold for features for
#' var_p_resid method
#' @param var_number (var_p_resid) number of top variance features for
#' var_p_resid method
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
        var_threshold = 1.5,
        var_number = NULL,
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

    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
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
        readGiottoInstructions(gobject, param = "show_plot"),
        show_plot
    )
    save_plot <- ifelse(is.na(save_plot),
        readGiottoInstructions(gobject, param = "save_plot"),
        save_plot
    )
    return_plot <- ifelse(is.na(return_plot),
        readGiottoInstructions(gobject, param = "return_plot"),
        return_plot
    )


    # method to use
    method <- match.arg(
        method,
        choices = c("cov_groups", "cov_loess", "var_p_resid")
    )

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
                analyzeParam("var", use_parallel = use_parallel)
            )
            if (!is.null(var_number) && is.numeric(var_number)) {
                dt[, selected := seq_len(.N)]
                dt[, selected := ifelse(selected <= var_number, "yes", "no")]
            } else {
                dt[, selected := ifelse(var >= var_threshold, "yes", "no")]
            }
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
            dt[, selected := ifelse(cov_group_zscore > zscore_threshold,
                                     "yes", "no")]
            dt
        },
        "cov_loess" = {
            dt <- analyzeData(
                expr_values,
                analyzeParam("cov_loess",
                    detection_threshold = expression_threshold,
                    use_parallel = use_parallel)
            )
            dt[, selected := ifelse(cov_diff > difference_in_cov, "yes", "no")]
            dt
        }
    )

    # Plot generation. The analyzeData methods drop intermediate columns
    # (expr_groups, pred_cov_feats, rank) used only for plotting; re-add
    # them locally so the existing plot helpers can be used unchanged.
    pl <- NULL
    if (want_plot) {
        if (method == "var_p_resid") {
            feat_in_cells_detected[, rank := seq_len(.N)]
            pl <- .create_calc_var_hvf_plot(feat_in_cells_detected)
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


.create_calc_var_hvf_plot <- function(dt_res) {
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_point(
        data = dt_res, 
        GiottoVisuals::aes_string2(x = "rank", y = "var", color = "selected")
    )
    pl <- pl + ggplot2::scale_x_reverse()
    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(
        axis.title = ggplot2::element_text(size = 14),
        axis.text = ggplot2::element_text(size = 12)
    )
    pl <- pl + ggplot2::scale_color_manual(
        values = c(no = "lightgrey", yes = "orange"),
        guide = ggplot2::guide_legend(
            title = "HVF",
            override.aes = list(size = 5)
        )
    )
    pl <- pl + ggplot2::labs(x = "feature rank", y = "variance")
    pl
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
setMethod("analyzeData",
    signature(x = "allMatrix", param = "featStatsParam"),
    function(x, param, ...) {
        mean_expr_det <- NULL
        det_thresh <- param$detection_threshold
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

        calc_fun <- if (isTRUE(param$use_parallel)) {
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

        calc_fun <- if (isTRUE(param$use_parallel)) {
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
        if (inherits(x, "IterableMatrix")) {
            scores <- sort(BPCells::rowVars(x), decreasing = TRUE)
            return(data.table::data.table(feats = names(scores), var = scores))
        }
        if (isTRUE(param$use_parallel)) {
            scores <- future.apply::future_apply(
                X = x, MARGIN = 1, FUN = var, future.seed = TRUE
            )
        } else {
            scores <- apply(X = x, MARGIN = 1, FUN = var)
        }
        scores <- sort(scores, decreasing = TRUE)
        data.table::data.table(feats = names(scores), var = scores)
    }
)

