.calc_cov_group_hvf <- function(
        feat_in_cells_detected,
        nr_expression_groups = 20,
        zscore_threshold = 1,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL) {
    # NSE vars
    cov_group_zscore <- cov <- selected <- mean_expr <- NULL

    steps <- 1 / nr_expression_groups
    prob_sequence <- seq(0, 1, steps)
    prob_sequence[length(prob_sequence)] <- 1
    expr_group_breaks <- stats::quantile(
        feat_in_cells_detected$mean_expr,
        probs = prob_sequence
    )

    ## remove zero's from cuts if there are too many and make first group zero
    if (any(duplicated(expr_group_breaks))) {
        m_expr_vector <- feat_in_cells_detected$mean_expr
        expr_group_breaks <- stats::quantile(
            m_expr_vector[m_expr_vector > 0],
            probs = prob_sequence
        )
        expr_group_breaks[[1]] <- 0
    }

    expr_groups <- cut(
        x = feat_in_cells_detected$mean_expr,
        breaks = expr_group_breaks,
        labels = paste0("group_", seq_len(nr_expression_groups)),
        include.lowest = TRUE
    )
    feat_in_cells_detected[, expr_groups := expr_groups]
    feat_in_cells_detected[, cov_group_zscore := scale(cov), by = expr_groups]
    feat_in_cells_detected[, selected := ifelse(
        cov_group_zscore > zscore_threshold, "yes", "no"
    )]

    if (any(isTRUE(show_plot), isTRUE(return_plot), isTRUE(save_plot))) {
        pl <- .create_cov_group_hvf_plot(
            feat_in_cells_detected, nr_expression_groups
        )

        return(list(dt = feat_in_cells_detected, pl = pl))
    } else {
        return(list(dt = feat_in_cells_detected))
    }
}







.calc_cov_loess_hvf <- function(
        feat_in_cells_detected,
        difference_in_cov = 0.1,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL) {
    # NSE vars
    cov_diff <- pred_cov_feats <- selected <- NULL

    # create loess regression
    loess_formula <- paste0("cov~log(mean_expr)")
    var_col <- "cov"

    loess_model_sample <- stats::loess(
        loess_formula,
        data = feat_in_cells_detected
    )
    feat_in_cells_detected$pred_cov_feats <- stats::predict(
        loess_model_sample,
        newdata = feat_in_cells_detected
    )
    feat_in_cells_detected[, cov_diff := get(var_col) - pred_cov_feats,
        by = seq_len(nrow(feat_in_cells_detected))
    ]
    data.table::setorder(feat_in_cells_detected, -cov_diff)
    feat_in_cells_detected[, selected := ifelse(
        cov_diff > difference_in_cov, "yes", "no"
    )]

    if (any(isTRUE(show_plot), isTRUE(return_plot), isTRUE(save_plot))) {
        pl <- .create_cov_loess_hvf_plot(
            feat_in_cells_detected, difference_in_cov, var_col
        )

        return(list(dt = feat_in_cells_detected, pl = pl))
    } else {
        return(list(dt = feat_in_cells_detected))
    }
}



.calc_var_hvf <- function(
        scaled_matrix,
        var_threshold = 1.5,
        var_number = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        use_parallel = FALSE) {
    # NSE vars
    var <- selected <- NULL

    if (isTRUE(use_parallel)) {
        test <- apply(X = scaled_matrix, MARGIN = 1, FUN = function(x) var(x))
    } else {
        test <- future.apply::future_apply(
            X = scaled_matrix, MARGIN = 1, FUN = function(x) var(x),
            future.seed = TRUE
        )
    }

    test <- sort(test, decreasing = TRUE)

    dt_res <- data.table::data.table(feats = names(test), var = test)

    if (!is.null(var_number) & is.numeric(var_number)) {
        dt_res[, selected := seq_len(.N)]
        dt_res[, selected := ifelse(selected <= var_number, "yes", "no")]
    } else {
        dt_res[, selected := ifelse(var >= var_threshold, "yes", "no")]
    }


    if (isTRUE(show_plot) ||
        isTRUE(return_plot) ||
        isTRUE(save_plot)) {
        dt_res[, rank := seq_len(.N)]
        pl <- .create_calc_var_hvf_plot(dt_res)


        dt_res_final <- data.table::copy(dt_res)
        dt_res_final[, rank := NULL]

        return(list(dt = dt_res_final, pl = pl))
    } else {
        return(list(dt = dt_res))
    }
}


.calc_expr_general_stats <- function(expr_values, expression_threshold) {
    # NSE vars
    gini <- NULL

    ## create data.table with relevant statistics ##
    feat_in_cells_detected <- data.table::data.table(
        feats = rownames(expr_values),
        nr_cells = rowSums_flex(expr_values > expression_threshold),
        total_expr = rowSums_flex(expr_values),
        mean_expr = rowMeans_flex(expr_values),
        sd = unlist(apply(expr_values, 1, sd))
    )

    # calculate gini rowwise
    gini_level <- unlist(apply(expr_values, MARGIN = 1, mygini_fun))
    feat_in_cells_detected[, gini := gini_level]

    return(feat_in_cells_detected)
}


.calc_expr_cov_stats <- function(expr_values, expression_threshold) {
    # NSE vars
    cov <- sd <- mean_expr <- NULL

    # get general expression statistics and gini data.table
    feat_in_cells_detected <- .calc_expr_general_stats(
        expr_values, expression_threshold
    )

    # calculate cov using sd and mean_expr from general stats DT
    feat_in_cells_detected[, cov := (sd / mean_expr)]

    return(feat_in_cells_detected)
}


.calc_expr_cov_stats_parallel <- function(expr_values,
    expression_threshold,
    cores = GiottoUtils::determine_cores()) {
    # NSE vars
    cov <- sd <- mean_expr <- NULL

    # setup chunk rows to use for each parallel based on number of cores
    chunk_rows <- seq(nrow(expr_values)) %>%
        split(., cut(., cores))

    # params to pass into the future_lapply
    fparams <- list(
        calc_fun = .calc_expr_general_stats,
        expression_threshold = expression_threshold
    )

    # parallelized calculation of general stats
    chunk_stats_dt_list <- lapply_flex(
        chunk_rows,
        function(r_idx, fparams) {
            fparams$calc_fun(
                expr_values = expr_values[r_idx, ],
                expression_threshold = fparams$expression_threshold
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
        if (isTRUE(set_seed)) set.seed(seed = seed_number)

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
    # select function to use based on whether future parallelization is planned
    calc_cov_fun <- ifelse(
        use_parallel,
        .calc_expr_cov_stats_parallel,
        .calc_expr_cov_stats
    )

    results <- switch(method,
        "var_p_resid" = {
            .calc_var_hvf(
                scaled_matrix = expr_values,
                var_threshold = var_threshold,
                var_number = var_number,
                show_plot = show_plot,
                return_plot = return_plot,
                save_plot = save_plot,
                use_parallel = use_parallel
            )
        },
        "cov_groups" = {
            calc_cov_fun(expr_values, expression_threshold) %>%
                .calc_cov_group_hvf(
                    nr_expression_groups = nr_expression_groups,
                    zscore_threshold = zscore_threshold,
                    show_plot = show_plot,
                    return_plot = return_plot,
                    save_plot = save_plot
                )
        },
        "cov_loess" = {
            calc_cov_fun(expr_values, expression_threshold) %>%
                .calc_cov_loess_hvf(
                    difference_in_cov = difference_in_cov,
                    show_plot = show_plot,
                    return_plot = return_plot,
                    save_plot = save_plot
                )
        }
    )

    ## unpack results
    feat_in_cells_detected <- results[["dt"]]
    pl <- results[["pl"]]




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
        ggplot2::aes_string(x = "mean_expr", y = "cov", color = "selected")
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
        ggplot2::aes_string(
            x = "log(mean_expr)", y = var_col,
            color = "selected"
        )
    )
    pl <- pl + ggplot2::geom_line(
        data = feat_in_cells_detected,
        ggplot2::aes_string(x = "log(mean_expr)", y = "pred_cov_feats"),
        color = "blue"
    )
    hvg_line <- paste0("pred_cov_feats+", difference_in_cov)
    pl <- pl + ggplot2::geom_line(
        data = feat_in_cells_detected,
        ggplot2::aes_string(x = "log(mean_expr)", y = hvg_line), linetype = 2
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
        data = dt_res, aes_string(x = "rank", y = "var", color = "selected")
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
