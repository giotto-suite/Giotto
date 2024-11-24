### Filter Values ####




#' @title filterDistributions
#' @name filterDistributions
#' @description show gene or cell distribution after filtering on expression
#' threshold
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param method method to create distribution (see details)
#' @param expression_threshold threshold to consider a gene expressed
#' @param detection consider features (e.g. genes) or cells
#' @param plot_type type of plot
#' @param scale_y scale y-axis (e.g. "log"), NULL = no scaling
#' @param nr_bins number of bins for histogram plot
#' @param fill_color fill color for plots
#' @param scale_axis ggplot transformation for axis (e.g. log2)
#' @param axis_offset offset to be used together with the scaling
#' transformation
#' @param show_plot logical. show plot
#' @param return_plot logical. return ggplot object
#' @param save_plot logical. directly save the plot
#' @param save_param list of saving parameters from
#' [GiottoVisuals::all_plots_save_function]
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @returns ggplot object
#' @details
#' There are 3 ways to create a distribution profile and summarize it for
#' either the features or the cells (spatial units) \cr
#' \itemize{
#'   \item{1. threshold: calculate features that cross a thresold (default)}
#'   \item{2. sum: summarize the features, i.e. total of a feature}
#'   \item{3. mean: calculate mean of the features, i.e. average expression}
#' }
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' filterDistributions(g)
#' @export
filterDistributions <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = c("raw", "normalized", "scaled", "custom"),
    method = c("threshold", "sum", "mean"),
    expression_threshold = 1,
    detection = c("feats", "cells"),
    plot_type = c("histogram", "violin"),
    scale_y = NULL,
    nr_bins = 30,
    fill_color = "lightblue",
    scale_axis = "identity",
    axis_offset = 0,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "filterDistributions") {
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
        unique(c("raw", "normalized", "scaled", "custom", expression_values))
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )

    # plot distribution for feats or cells
    detection <- match.arg(detection, c("feats", "cells"))

    # method to calculate distribution
    method <- match.arg(method, c("threshold", "sum", "mean"))

    # plot type
    plot_type <- match.arg(plot_type, c("histogram", "violin"))

    # variables
    V1 <- NULL

    # for genes
    if (detection == "feats") {
        if (method == "threshold") {
            feat_detection_levels <- data.table::as.data.table(
                rowSums_flex(expr_values >= expression_threshold)
            )
            mytitle <- "feat detected in # of cells"
        } else if (method == "sum") {
            feat_detection_levels <- data.table::as.data.table(
                rowSums_flex(expr_values)
            )
            mytitle <- "total sum of feature detected in all cells"
        } else if (method == "mean") {
            feat_detection_levels <- data.table::as.data.table(
                rowMeans_flex(expr_values)
            )
            mytitle <- "average of feature detected in all cells"
        }

        y_title <- "count"
        if (!is.null(scale_y)) {
            feat_detection_levels[, V1 := do.call(what = scale_y, list(V1))]
            y_title <- paste(scale_y, y_title)
        }



        if (plot_type == "violin") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::geom_violin(
                data = feat_detection_levels,
                ggplot2::aes(x = "feats", y = V1 + axis_offset),
                fill = fill_color
            )
            pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
            pl <- pl + ggplot2::labs(y = mytitle, x = "")
        } else if (plot_type == "histogram") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::geom_histogram(
                data = feat_detection_levels,
                ggplot2::aes(x = V1 + axis_offset),
                color = "white", bins = nr_bins, fill = fill_color
            )
            pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
            pl <- pl + ggplot2::labs(x = mytitle, y = y_title)
        }

        # for cells
    } else if (detection == "cells") {
        if (method == "threshold") {
            cell_detection_levels <- data.table::as.data.table(
                colSums_flex(expr_values >= expression_threshold)
            )
            mytitle <- "feats detected per cell"
        } else if (method == "sum") {
            cell_detection_levels <- data.table::as.data.table(
                colSums_flex(expr_values)
            )
            mytitle <- "total features per cell"
        } else if (method == "mean") {
            cell_detection_levels <- data.table::as.data.table(
                colMeans_flex(expr_values)
            )
            mytitle <- "average number of features per cell"
        }

        y_title <- "count"
        if (!is.null(scale_y)) {
            cell_detection_levels[, V1 := do.call(what = scale_y, list(V1))]
            y_title <- paste(scale_y, y_title)
        }



        if (plot_type == "violin") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::geom_violin(
                data = cell_detection_levels,
                ggplot2::aes(x = "cells", y = V1 + axis_offset),
                fill = fill_color
            )
            pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
            pl <- pl + ggplot2::labs(y = mytitle, x = "")
        } else if (plot_type == "histogram") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::geom_histogram(
                data = cell_detection_levels,
                ggplot2::aes(x = V1 + axis_offset),
                color = "white", bins = nr_bins, fill = fill_color
            )
            pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
            pl <- pl + ggplot2::labs(x = mytitle, y = y_title)
        }
    }

    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}



#' @title filterCombinations
#' @name filterCombinations
#' @description Shows how many genes and cells are lost with combinations of
#' thresholds.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param expression_values expression values to use
#' @param expression_thresholds all thresholds to consider a gene expressed
#' @param feat_det_in_min_cells minimum # of cells that need to express a
#' feature
#' @param min_det_feats_per_cell minimum # of features that need to be
#' detected in a cell
#' @param scale_x_axis ggplot transformation for x-axis (e.g. log2)
#' @param x_axis_offset x-axis offset to be used together with the scaling
#' transformation
#' @param scale_y_axis ggplot transformation for y-axis (e.g. log2)
#' @param y_axis_offset y-axis offset to be used together with the scaling
#' transformation
#' @returns list of data.table and ggplot object
#' @details Creates a scatterplot that visualizes the number of genes and
#' cells that are lost with a specific combination of a gene and cell
#' threshold given an arbitrary cutoff to call a gene expressed. This function
#' can be used to make an informed decision at the filtering step with
#' filterGiotto.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' filterCombinations(g)
#' @export
filterCombinations <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = c("raw", "normalized", "scaled", "custom"),
    expression_thresholds = c(1, 2),
    feat_det_in_min_cells = c(5, 50),
    min_det_feats_per_cell = c(200, 400),
    scale_x_axis = "identity",
    x_axis_offset = 0,
    scale_y_axis = "identity",
    y_axis_offset = 0,
    show_plot = TRUE,
    return_plot = FALSE,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "filterCombinations") {
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
        unique(c("raw", "normalized", "scaled", "custom", expression_values))
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values
    )[]

    # feat and cell minimums need to have the same length
    if (length(feat_det_in_min_cells) != length(min_det_feats_per_cell)) {
        stop("\n feat_det_in_min_cells and min_det_feats_per_cell need to be
            the same size \n")
    }

    # compute the number of removed feats and cells
    result_list <- list()
    for (thresh_i in seq_along(expression_thresholds)) {
        threshold <- expression_thresholds[thresh_i]

        det_feats_res <- list()
        det_cells_res <- list()
        for (combn_i in seq_along(feat_det_in_min_cells)) {
            min_cells_for_feat <- feat_det_in_min_cells[combn_i]
            min_feats_per_cell <- min_det_feats_per_cell[combn_i]


            # first remove feats
            filter_index_feats <- rowSums_flex(
                expr_values >= threshold
            ) >= min_cells_for_feat
            removed_feats <- length(filter_index_feats[
                filter_index_feats == FALSE
            ])
            det_cells_res[[combn_i]] <- removed_feats

            # then remove cells
            filter_index_cells <- colSums_flex(expr_values[
                filter_index_feats,
            ] >= threshold) >= min_feats_per_cell
            removed_cells <- length(filter_index_cells[
                filter_index_cells == FALSE
            ])
            det_feats_res[[combn_i]] <- removed_cells
        }

        temp_dt <- data.table::data.table(
            "threshold" = threshold,
            removed_feats = unlist(det_cells_res),
            removed_cells = unlist(det_feats_res)
        )

        result_list[[thresh_i]] <- temp_dt
    }

    result_DT <- do.call("rbind", result_list)

    # data.table variables
    feat_detected_in_min_cells <- min_detected_feats_per_cell <-
        combination <- NULL

    result_DT[["feat_detected_in_min_cells"]] <- feat_det_in_min_cells
    result_DT[["min_detected_feats_per_cell"]] <- min_det_feats_per_cell
    result_DT[["combination"]] <- paste0(
        result_DT$feat_detected_in_min_cells, "-",
        result_DT$min_detected_feats_per_cell
    )

    result_DT <- result_DT[, .(
        threshold,
        feat_detected_in_min_cells,
        min_detected_feats_per_cell,
        combination,
        removed_feats,
        removed_cells
    )]

    maximum_x_value <- max(result_DT[["removed_cells"]], na.rm = TRUE)
    maximum_y_value <- max(result_DT[["removed_feats"]], na.rm = TRUE)

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()
    pl <- pl + ggplot2::geom_line(data = result_DT, aes(
        x = removed_cells + x_axis_offset,
        y = removed_feats + y_axis_offset,
        group = as.factor(threshold)
    ), linetype = 2)
    pl <- pl + ggplot2::geom_point(data = result_DT, aes(
        x = removed_cells + x_axis_offset,
        y = removed_feats + y_axis_offset,
        color = as.factor(threshold)
    ))
    pl <- pl + scale_color_discrete(
        guide = guide_legend(title = "threshold(s)")
    )
    pl <- pl + geom_text_repel(data = result_DT, aes(
        x = removed_cells + x_axis_offset,
        y = removed_feats + y_axis_offset,
        label = combination
    ))
    pl <- pl + ggplot2::scale_x_continuous(
        trans = scale_x_axis, limits = c(0, maximum_x_value)
    )
    pl <- pl + ggplot2::scale_y_continuous(
        trans = scale_y_axis, limits = c(0, maximum_y_value)
    )
    pl <- pl + ggplot2::labs(
        x = "number of removed cells", y = "number of removed feats"
    )


    return(plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = list(results = result_DT, ggplot = pl)
    ))
}


#' @title filterGiotto
#' @name filterGiotto
#' @description filter Giotto object based on expression threshold
#' @param gobject giotto object
#' @param spat_unit character. spatial unit. If more than one is provided then
#' the first will be filtered, the filtering results will be applied across the
#' other spat_units provided
#' @param feat_type character. feature type. If more than one is provided then
#' the first will be filtered, the filtering results will be applied across the
#' other feat_types provided.
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param feat_det_in_min_cells minimum # of cells that need to express a
#' feature
#' @param min_det_feats_per_cell minimum # of features that need to be detected
#' in a cell
#' @param all_spat_units deprecated. Use spat_unit_fsub = ":all:"
#' @param all_feat_types deprecated. Use feat_type_ssub = ":all:"
#' @param spat_unit_fsub character vector. (default = ':all:') limit features
#' to remove results to selected spat_units
#' @param feat_type_ssub character vector. (default = ':all:') limit cells to
#' remove results to selected feat_types
#' @param poly_info polygon information to use
#' @param tag_cells tag filtered cells in metadata vs. remove cells
#' @param tag_cell_name column name for tagged cells in metadata
#' @param tag_feats tag features in metadata vs. remove features
#' @param tag_feats_name column name for tagged features in metadata
#' @param verbose verbose
#'
#' @returns giotto object
#' @details The function \code{\link{filterCombinations}} can be used to
#' explore the effect of different parameter values.
#' Please note that this function filters data in a predefined order, features,
#' then cells.
#' After filtering in this order, certain features may be left over in the
#' metadata with a corresponding number of cells which is less than that of
#' the threshold value of cells,
#' feat_det_in_min_cells. This behavior is explained in detail here:
#' \url{https://github.com/drieslab/Giotto/issues/500#issuecomment-1396083446}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' filterGiotto(g)
#' @export
filterGiotto <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("raw", "normalized", "scaled", "custom"),
    expression_threshold = 1,
    feat_det_in_min_cells = 100,
    min_det_feats_per_cell = 100,
    spat_unit_fsub = ":all:",
    feat_type_ssub = ":all:",
    all_spat_units = NULL,
    all_feat_types = NULL,
    poly_info = NULL,
    tag_cells = FALSE,
    tag_cell_name = "tag",
    tag_feats = FALSE,
    tag_feats_name = "tag",
    verbose = TRUE) {
    # data.table vars
    cell_ID <- feat_ID <- NULL

    # handle deprecations
    if (!is.null(all_spat_units)) {
        if (all_spat_units) {
            spat_unit_fsub <- ":all:"
        } else {
            spat_unit_fsub <- spat_unit
        }

        warning(wrap_txt(
            'filterGiotto:
      all_spat_units param is deprecated.
      Please use spat_unit_fsub = \":all:\" instead. (this is the default)'
        ))
    }
    if (!is.null(all_feat_types)) {
        if (all_feat_types) {
            feat_type_ssub <- ":all:"
        } else {
            feat_type_ssub <- feat_type
        }

        warning(wrap_txt(
            'filterGiotto: all_feat_types param is deprecated.
            Please use feat_type_ssub = \":all:\" instead.
            (this is the default)'
        ))
    }


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
    # set poly_info
    if (is.null(poly_info)) {
        poly_info <- spat_unit
    }

    if (verbose && length(spat_unit) > 1L) {
        wrap_msg(
            "More than one spat_unit provided.\n",
            paste0("[", spat_unit[[1L]], "]"),
            "filtering results will be applied across spat_units:", spat_unit
        )
    }
    if (verbose && length(feat_type) > 1L) {
        wrap_msg(
            "More than one feat_type provided.\n",
            paste0("[", feat_type[[1L]], "]"),
            "filtering results will be applied across spat_units:", feat_type
        )
    }


    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("raw", "normalized", "scaled", "custom", expression_values))
    )

    # get expression values to perform filtering on
    # Only the first spat_unit and feat_type provided are filtered.
    # IF there are additional spat_units and feat_types provided, then the
    # filtering
    # results from this round will be applied to the other provided spat_units
    # and feat_types as well.
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit[[1L]],
        feat_type = feat_type[[1L]],
        values = values,
        output = "matrix"
    )

    # approach:
    # 1. first remove genes that are not frequently detected
    # 2. then remove cells that do not have sufficient detected genes

    ## filter features
    filter_index_feats <- rowSums_flex(
        expr_values >= expression_threshold
    ) >= feat_det_in_min_cells
    selected_feat_ids <- names(filter_index_feats[filter_index_feats == TRUE])



    ## filter cells
    filter_index_cells <- colSums_flex(expr_values[
        filter_index_feats,
    ] >= expression_threshold) >= min_det_feats_per_cell
    selected_cell_ids <- names(filter_index_cells[filter_index_cells == TRUE])



    # update cell metadata
    if (isTRUE(tag_cells)) {
        cell_meta <- getCellMetadata(gobject = gobject, copy_obj = TRUE)
        cell_meta[][, c(tag_cell_name) := ifelse(
            cell_ID %in% selected_cell_ids, 0, 1
        )]
        gobject <- setCellMetadata(
            gobject = gobject, x = cell_meta, initialize = FALSE
        )

        # set selected cells back to all cells
        selected_cell_ids <- names(filter_index_cells)
    }

    if (isTRUE(tag_feats)) {
        feat_meta <- getFeatureMetadata(gobject = gobject, copy_obj = TRUE)
        feat_meta[][, c(tag_feats_name) := ifelse(
            feat_ID %in% selected_feat_ids, 0, 1
        )]
        gobject <- setFeatureMetadata(
            gobject = gobject, x = feat_meta, initialize = FALSE
        )

        # set selected feats back to all feats
        selected_feat_ids <- names(filter_index_feats)
    }



    # update feature metadata
    newGiottoObject <- subsetGiotto(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        cell_ids = selected_cell_ids,
        feat_ids = selected_feat_ids,
        spat_unit_fsub = spat_unit_fsub,
        feat_type_ssub = feat_type_ssub,
        poly_info = poly_info,
        verbose = verbose
    )

    ## print output ##
    removed_feats <- length(filter_index_feats[filter_index_feats == FALSE])
    total_feats <- length(filter_index_feats)

    removed_cells <- length(filter_index_cells[filter_index_cells == FALSE])
    total_cells <- length(filter_index_cells)

    if (isTRUE(verbose)) {
        cat("\n")
        cat("Feature type: ", feat_type, "\n")

        if (isTRUE(tag_cells)) {
            cat(
                "Number of cells tagged: ", removed_cells, " out of ",
                total_cells, "\n"
            )
        } else {
            cat(
                "Number of cells removed: ", removed_cells, " out of ",
                total_cells, "\n"
            )
        }

        if (isTRUE(tag_feats)) {
            cat(
                "Number of feats tagged: ", removed_feats, " out of ",
                total_feats, "\n"
            )
        } else {
            cat(
                "Number of feats removed: ", removed_feats, " out of ",
                total_feats, "\n"
            )
        }
    }


    ## update parameters used ##

    # Do not update downstream of processGiotto
    # Parameters will be updated within processGiotto
    try(
        {
            upstream_func <- sys.call(-2)
            fname <- as.character(upstream_func[[1]])
            if (fname == "processGiotto") {
                return(newGiottoObject)
            }
        },
        silent = TRUE
    )


    # If this function call is not downstream of processGiotto, update normally
    newGiottoObject <- update_giotto_params(
        newGiottoObject,
        description = "_filter"
    )

    return(newGiottoObject)
}
