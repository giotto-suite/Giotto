## Giotto auxiliary functions ####






#' @title Adjust expression values
#' @name adjustGiottoMatrix
#' @description Adjust expression values to account for known batch effects or
#' technological covariates.
#' @inheritParams data_access_params
#' @param expression_values expression values to use
#' @param batch_columns metadata columns that represent different
#' batch (max = 2)
#' @param covariate_columns metadata columns that represent covariates to
#' regress out
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param update_slot expression slot that will be updated (default = custom)
#' @returns giotto object or exprObj
#' @details This function implements the \code{\link[limma]{removeBatchEffect}}
#' function to remove known batch effects and to adjust expression values
#' according to provided covariates.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' adjustGiottoMatrix(g, covariate_columns = "leiden_clus")
#' @export
adjustGiottoMatrix <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    batch_columns = NULL,
    covariate_columns = NULL,
    return_gobject = TRUE,
    update_slot = c("custom")) {
    # Catch for both batch and covariate being null
    if (is.null(batch_columns) & is.null(covariate_columns)) {
        stop("Metadata for either different batches or covariates must be
            provided.")
    }

    package_check("limma")

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

    # metadata
    cell_metadata <- getCellMetadata(
        gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        output = "data.table",
        copy_obj = TRUE
    )

    if (!is.null(batch_columns)) {
        if (!all(batch_columns %in% colnames(cell_metadata))) {
            stop("batch column name(s) were not found in the cell metadata")
        }
    }

    if (!is.null(covariate_columns)) {
        if (!all(covariate_columns %in% colnames(cell_metadata))) {
            stop("covariate column name(s) were not found in the cell metadata")
        }
    }

    update_slot <- match.arg(
        update_slot, c("normalized", "scaled", "custom", update_slot)
    )

    # expression values to be used
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "exprObj"
    )


    # batch columns
    if (!is.null(batch_columns)) {
        batch_column_1 <- cell_metadata[[batch_columns[1]]]
        if (length(batch_columns) > 1) {
            batch_column_2 <- cell_metadata[[batch_columns[2]]]
        } else {
            batch_column_2 <- NULL
        }
    } else {
        batch_column_1 <- NULL
        batch_column_2 <- NULL
    }

    # covariate columns
    if (!is.null(covariate_columns)) {
        covariates <- as.matrix(
            cell_metadata[, covariate_columns, with = FALSE]
        )
    } else {
        covariates <- NULL
    }



    # TODO: implement ResidualMatrix to work with a delayed matrix
    adjusted_matrix <- limma::removeBatchEffect(
        x = expr_data[],
        batch = batch_column_1,
        batch2 = batch_column_2,
        covariates = covariates
    )

    if (return_gobject == TRUE) {
        adjusted_matrix <- create_expr_obj(
            name = update_slot,
            exprMat = adjusted_matrix,
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = expr_data@provenance
        )

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- setGiotto(gobject, adjusted_matrix)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        ## update parameters used ##

        # Do not update downstream of processGiotto
        # Parameters will be updated within processGiotto
        try(
            {
                test <- sys.call(-2)
                fname <- as.character(test[[1]])
                if (fname == "processGiotto") {
                    return(gobject)
                }
            },
            silent = TRUE
        )


        # If this function call is not downstream of processGiotto, update
        # normally
        gobject <- update_giotto_params(gobject, description = "_adj_matrix")

        return(gobject)
    } else {
        return(adjusted_matrix)
    }
}


#' @title processGiotto
#' @name processGiotto
#' @description Wrapper for the different Giotto object processing functions
#' @param gobject giotto object
#' @param filter_params additional parameters to filterGiotto
#' @param norm_params additional parameters to normalizeGiotto
#' @param stat_params additional parameters to addStatistics
#' @param adjust_params additional parameters to adjustGiottoMatrix; set to
#' NULL if not required
#' @param verbose be verbose (default is TRUE)
#' @returns giotto object
#' @details See \code{\link{filterGiotto}}, \code{\link{normalizeGiotto}},
#' \code{\link{addStatistics}}, and \code{\link{adjustGiottoMatrix}}. For more
#' information about the different parameters in each step. If you do not
#' provide them it will use the default values. If no adjustment is required,
#' adjust_params must be set to NULL
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' processGiotto(
#'     gobject = g,
#'     adjust_params = list(covariate_columns = "leiden_clus")
#' )
#' @export
processGiotto <- function(gobject,
    filter_params = list(),
    norm_params = list(),
    stat_params = list(),
    adjust_params = list(),
    verbose = TRUE) {
    # filter Giotto
    vmsg(.v = verbose, "1. start filter step")
    if (!inherits(filter_params, "list")) {
        stop("filter_params need to be a list of parameters for filterGiotto")
    }
    gobject <- do.call("filterGiotto", c(gobject = gobject, filter_params))

    # normalize Giotto
    vmsg(.v = verbose, "2. start normalization step")
    if (!inherits(norm_params, "list")) {
        stop("norm_params need to be a list of parameters for normalizeGiotto")
    }
    gobject <- do.call("normalizeGiotto", c(gobject = gobject, norm_params))

    # add Statistics
    if (verbose == TRUE) message("3. start cell and gene statistics step")
    if (!inherits(stat_params, "list")) {
        stop("stat_params need to be a list of parameters for addStatistics ")
    }
    stat_params[["return_gobject"]] <- TRUE # force this to be true
    gobject <- do.call("addStatistics", c(gobject = gobject, stat_params))

    # adjust Giotto, if applicable
    if (length(adjust_params) > 0L) {
        if (verbose == TRUE) message("4. start adjusted matrix step")
        if (!inherits(adjust_params, "list")) {
            stop("adjust_params need to be a list of parameters for
                adjustGiottoMatrix")
        }
        adjust_params[["return_gobject"]] <- TRUE # force this to be true
        gobject <- do.call(
            "adjustGiottoMatrix", c(gobject = gobject, adjust_params)
        )
    }

    gobject <- update_giotto_params(gobject, description = "_process")

    return(gobject)
}










## * ####
## Feature & Cell statistics ####





#' @title Add feature statistics
#' @name addFeatStatistics
#' @description Adds feature statistics to the giotto object
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose be verbose
#' @returns giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to feature metadata:
#' \itemize{
#'   \item{\strong{nr_cells:} Denotes in how many cells the feature is
#'   detected}
#'   \item{\strong{per_cells:} Denotes in what percentage of cells the feature
#'   is detected}
#'   \item{\strong{total_expr:} Shows the total sum of feature expression in
#'   all cells}
#'   \item{\strong{mean_expr:} Average feature expression in all cells}
#'   \item{\strong{mean_expr_det:} Average feature expression in cells with
#'   detectable levels of the gene}
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' addFeatStatistics(g)
#' @export
addFeatStatistics <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    detection_threshold = 0,
    return_gobject = TRUE,
    verbose = TRUE) {
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
    expression_values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = expression_values,
        output = "exprObj",
        set_defaults = FALSE
    )

    # calculate stats
    feat_stats <- data.table::data.table(
        feats = rownames(expr_data[]),
        nr_cells = rowSums_flex(expr_data[] > detection_threshold),
        perc_cells = (rowSums_flex(expr_data[] > detection_threshold) /
            ncol(expr_data[])) * 100,
        total_expr = rowSums_flex(expr_data[]),
        mean_expr = rowMeans_flex(expr_data[])
    )

    # data.table variables
    mean_expr_det <- NULL

    mean_expr_detected <- .mean_expr_det_test(
        expr_data[],
        detection_threshold = detection_threshold
    )
    feat_stats[, mean_expr_det := mean_expr_detected]


    if (return_gobject == TRUE) {
        # remove previous statistics
        feat_metadata <- getFeatureMetadata(
            gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "featMetaObj",
            copy_obj = TRUE,
            set_defaults = FALSE
        )

        if (isS4(expr_data)) {
            if (!identical(expr_data@provenance, feat_metadata@provenance)) {
                warning("expression and feature metadata provenance mismatch")
            }
        }


        metadata_names <- colnames(feat_metadata[])

        if ("nr_cells" %in% metadata_names) {
            vmsg(
                .v = verbose, "feat statistics has already been applied",
                "once; overwriting"
            )
            feat_metadata[][, c(
                "nr_cells", "perc_cells", "total_expr", "mean_expr",
                "mean_expr_det"
            ) := NULL]
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- setGiotto(gobject, feat_metadata, verbose = FALSE)
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }

        gobject <- addFeatMetadata(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            new_metadata = feat_stats,
            by_column = TRUE,
            column_feat_ID = "feats"
        )

        ## update parameters used ##

        # parent function name
        cl <- sys.call(-1)

        # Do not update downstream of processGiotto
        # Parameters will be updated within processGiotto
        try(
            {
                upstream_func <- sys.call(-3)
                fname <- as.character(upstream_func[[1]])
                if (fname == "processGiotto") {
                    return(gobject)
                }
            },
            silent = TRUE
        )


        # If this function call is not downstream of processGiotto, update
        # normally
        if (is.null(cl)) {
            gobject <- update_giotto_params(gobject,
                description = "_feat_stats"
            )
        } else {
            fname <- as.character(cl[[1]])
            if (fname == "addStatistics") {
                gobject <- update_giotto_params(gobject,
                    description = "_feat_stats",
                    toplevel = 3
                )
            } else {
                gobject <- update_giotto_params(gobject,
                    description = "_feat_stats"
                )
            }
        }


        return(gobject)
    } else {
        return(feat_stats)
    }
}





#' @title addCellStatistics
#' @name addCellStatistics
#' @description adds cells statistics to the giotto object
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose be verbose
#' @returns giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to cell metadata:
#' \itemize{
#'   \item{\strong{nr_feats:} Denotes in how many features are detected per
#'   cell}
#'   \item{\strong{perc_feats:} Denotes what percentage of features is
#'   detected per cell}
#'   \item{\strong{total_expr:} Shows the total sum of feature expression per
#'   cell}
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' addCellStatistics(g)
#' @export
addCellStatistics <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    detection_threshold = 0,
    return_gobject = TRUE,
    verbose = TRUE) {
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
    expression_values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = expression_values,
        output = "exprObj",
        set_defaults = FALSE
    )

    # calculate stats

    cell_stats <- data.table::data.table(
        cells = colnames(expr_data[]),
        nr_feats = colSums_flex(expr_data[] > detection_threshold),
        perc_feats = (colSums_flex(expr_data[] > detection_threshold) /
            nrow(expr_data[])) * 100,
        total_expr = colSums_flex(expr_data[])
    )

    if (return_gobject == TRUE) {
        # remove previous statistics
        cell_metadata <- getCellMetadata(
            gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "cellMetaObj",
            copy_obj = TRUE,
            set_defaults = FALSE
        )

        if (isS4(expr_data)) {
            if (!identical(expr_data@provenance, cell_metadata@provenance)) {
                warning("expression and feature metadata provenance mismatch")
            }
        }

        metadata_names <- colnames(cell_metadata[])
        if ("nr_feats" %in% metadata_names) {
            vmsg(
                .v = verbose, "cells statistics has already been applied",
                "once; overwriting"
            )
            cell_metadata[][, c("nr_feats", "perc_feats", "total_expr") := NULL]
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject <- setGiotto(gobject, cell_metadata, verbose = FALSE)
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }




        gobject <- addCellMetadata(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            new_metadata = cell_stats,
            by_column = TRUE,
            column_cell_ID = "cells"
        )

        ## update parameters used ##

        # parent function name
        cl <- sys.call(-1)

        # Do not update downstream of processGiotto
        # Parameters will be updated within processGiotto
        try(
            {
                upstream_func <- sys.call(-3)
                fname <- as.character(upstream_func[[1]])
                if (fname == "processGiotto") {
                    return(gobject)
                }
            },
            silent = TRUE
        )

        # If this function call is not downstream of processGiotto, update
        # normally
        if (is.null(cl)) {
            gobject <- update_giotto_params(gobject,
                description = "_cell_stats"
            )
        } else {
            fname <- as.character(cl[[1]])
            if (fname == "addStatistics") {
                gobject <- update_giotto_params(gobject,
                    description = "_cell_stats",
                    toplevel = 3
                )
            } else {
                gobject <- update_giotto_params(gobject,
                    description = "_cell_stats"
                )
            }
        }


        return(gobject)
    } else {
        return(cell_stats)
    }
}


#' @title addStatistics
#' @name addStatistics
#' @description Adds feature and cell statistics to the giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a feature detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose be verbose
#' @returns giotto object if return_gobject = TRUE, else a list with results
#' @details See \code{\link{addFeatStatistics}} and
#' \code{\link{addCellStatistics}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' addStatistics(g)
#' @export
addStatistics <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    detection_threshold = 0,
    return_gobject = TRUE,
    verbose = TRUE) {
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

    # get feats statistics
    feat_stats <- addFeatStatistics(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        expression_values = expression_values,
        detection_threshold = detection_threshold,
        return_gobject = return_gobject,
        verbose = verbose
    )

    if (return_gobject == TRUE) {
        gobject <- feat_stats
    }

    # get cell statistics
    cell_stats <- addCellStatistics(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        expression_values = expression_values,
        detection_threshold = detection_threshold,
        return_gobject = return_gobject,
        verbose = verbose
    )

    if (return_gobject == TRUE) {
        gobject <- cell_stats
        return(gobject)
    } else {
        return(feat_stats = feat_stats, cell_stats = cell_stats)
    }
}




#' @title addFeatsPerc
#' @name addFeatsPerc
#' @description Calculates the total percentage of (normalized) counts for a
#' subset of selected genes
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param feats vector of selected features
#' @param vector_name column name as seen in \code{\link{pDataDT}}
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @returns giotto object if \code{return_gobject = TRUE}, else a vector with %
#' results
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' addFeatsPerc(g, feats = c("Gm19935", "9630013A20Rik", "2900040C04Rik"))
#' @export
addFeatsPerc <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    feats = NULL,
    vector_name = "feat_perc",
    return_gobject = TRUE) {
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

    # tests
    if (is.null(feats)) {
        stop("You need to provide a vector of feat names")
    }

    if (!methods::is(gobject, "giotto")) {
        stop("You need to provide a giotto object")
    }


    # expression values to be used
    expression_values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = expression_values,
        output = "matrix"
    )


    totalsum <- colSums_flex(expr_data)
    feat_sum <- colSums_flex(expr_data[rownames(expr_data) %in% feats, ])
    perc_feats <- round((feat_sum / totalsum) * 100, 2)

    if (return_gobject) {
        temp_gobj <- addCellMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            new_metadata = perc_feats,
            vector_name = vector_name,
            by_column = TRUE
        )

        ## update parameters used ##
        temp_gobj <- update_giotto_params(temp_gobj,
            description = "_feats_perc"
        )

        return(temp_gobj)
    } else {
        return(perc_feats)
    }
}





## * ####
## Giotto auxiliary functions ####




#' @title Find network neighbors
#' @name findNetworkNeighbors
#' @description Find the spatial neighbors for a selected group of cells within
#' the selected spatial network.
#' @param gobject Giotto object
#' @param spat_unit spatial unit
#' @param spatial_network_name name of spatial network
#' @param source_cell_ids cell ids for which you want to know the spatial
#' neighbors
#' @param name name of the results
#' @returns data.table
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findNetworkNeighbors(
#'     gobject = g, spatial_network_name = "spatial_network",
#'     source_cell_ids = c("AACTCGATGGCGCAGT-1", "GGCTGGCTAGCTTAAA-1")
#' )
#' @export
findNetworkNeighbors <- function(gobject,
    spat_unit = NULL,
    spatial_network_name = NULL,
    source_cell_ids = NULL,
    name = "nb_cells") {
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )

    # get spatial network
    if (!is.null(spatial_network_name)) {
        spatial_network <- getSpatialNetwork(gobject,
            spat_unit = spat_unit,
            name = spatial_network_name,
            output = "networkDT"
        )
    } else {
        stop("You need to select a spatial network")
    }

    # source cell ids that are found back
    all_cell_ids <- gobject@cell_ID[[spat_unit]]
    source_cells <- all_cell_ids[all_cell_ids %in% source_cell_ids]

    if (length(source_cells) == 0) {
        stop("No source cell ids were selected or found")
    }


    full_network_DT <- convert_to_full_spatial_network(spatial_network)
    potential_target_cells <- full_network_DT[
        source %in% source_cells
    ][["target"]]
    source_and_target_cells <- potential_target_cells[
        potential_target_cells %in% source_cells
    ]
    target_cells <- potential_target_cells[
        !potential_target_cells %in% source_and_target_cells
    ]

    cell_meta <- pDataDT(gobject)

    # data.table variables
    nb_cells <- cell_ID <- NULL

    cell_meta[, nb_cells := ifelse(cell_ID %in% source_and_target_cells, "both",
        ifelse(cell_ID %in% source_cells, "source",
            ifelse(cell_ID %in% target_cells, "neighbor", "others")
        )
    )]
    nb_annot <- cell_meta[, c("cell_ID", "nb_cells"), with = FALSE]
    data.table::setnames(nb_annot, "nb_cells", name)

    return(nb_annot)
}


# internals ####



#' @title Mean expression detected test
#' @param mymatrix matrix of expression info
#' @param detection_threshold detection threshold. Defaults to 1 count.
#' @returns numeric
#' @keywords internal
#' @noRd
.mean_expr_det_test <- function(mymatrix, detection_threshold = 1) {
    unlist(apply(X = mymatrix, MARGIN = 1, FUN = function(x) {
        detected_x <- x[x > detection_threshold]
        mean(detected_x)
    }))
}
