#' @title findScranMarkers
#' @name findScranMarkers
#' @description Identify marker genes for all or selected clusters based on
#' scran's implementation of findMarkers.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values gene expression values to use
#' @param cluster_column clusters to use
#' @param subset_clusters selection of clusters to compare
#' @param group_1 group 1 cluster IDs from cluster_column for pairwise
#' comparison
#' @param group_1_name custom name for group_1 clusters
#' @param group_2 group 2 cluster IDs from cluster_column for pairwise
#' comparison
#' @param group_2_name custom name for group_2 clusters
#' @param verbose be verbose (default = FALSE)
#' @param ... additional parameters for the findMarkers function in scran
#' @returns data.table with marker genes
#' @details This is a minimal convenience wrapper around
#' the \code{\link[scran]{findMarkers}} function from the scran package.
#'
#' To perform differential expression between custom selected groups of cells
#' you need to specify the cell_ID column to parameter \emph{cluster_column}
#' and provide the individual cell IDs to the parameters \emph{group_1} and
#' \emph{group_2}
#'
#' By default group names will be created by pasting the different id names
#' within each selected group.
#' When you have many different ids in a single group
#' it is recommend to provide names for both groups to \emph{group_1_name} and
#' \emph{group_2_name}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findScranMarkers(g, cluster_column = "leiden_clus")
#' @export
findScranMarkers <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        cluster_column,
        subset_clusters = NULL,
        group_1 = NULL,
        group_1_name = NULL,
        group_2 = NULL,
        group_2_name = NULL,
        verbose = TRUE,
        ...) {
    # verify if optional package is installed
    package_check(pkg_name = "scran", repository = "Bioc")


    # print message with information #
    if (isTRUE(verbose)) {
        message("Using 'Scran' to detect marker genes. If used in published
        research, please cite:
        Lun ATL, McCarthy DJ, Marioni JC (2016).
        'A step-by-step workflow for low-level analysis of single-cell RNA-seq
        data with Bioconductor.'
        F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.")
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

    # expression data
    values <- match.arg(
        expression_values,
        choices = unique(c(
            "normalized", "scaled", "custom",
            expression_values
        ))
    )
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )

    # cluster column
    cell_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )
    if (!cluster_column %in% colnames(cell_metadata)) {
        stop("cluster column not found")
    }

    # subset expression data
    if (!is.null(subset_clusters)) {
        cell_metadata <- cell_metadata[get(cluster_column) %in% subset_clusters]
        subset_cell_IDs <- cell_metadata[["cell_ID"]]
        expr_data <- expr_data[, colnames(expr_data) %in% subset_cell_IDs]
    } else if (!is.null(group_1) & !is.null(group_2)) {
        cell_metadata <- cell_metadata[
            get(cluster_column) %in% c(group_1, group_2)
        ]

        # create new pairwise group
        if (!is.null(group_1_name)) {
            if (!is.character(group_1_name)) {
                stop("group_1_name needs to be a character")
            }
            group_1_name <- group_1_name
        } else {
            group_1_name <- paste0(group_1, collapse = "_")
        }

        if (!is.null(group_2_name)) {
            if (!is.character(group_2_name)) {
                stop("group_2_name needs to be a character")
            }
            group_2_name <- group_2_name
        } else {
            group_2_name <- paste0(group_2, collapse = "_")
        }


        # data.table variables
        pairwise_select_comp <- NULL

        cell_metadata[, pairwise_select_comp := ifelse(
            get(cluster_column) %in% group_1, group_1_name, group_2_name
        )]

        cluster_column <- "pairwise_select_comp"

        # expression data
        subset_cell_IDs <- cell_metadata[["cell_ID"]]
        expr_data <- expr_data[, colnames(expr_data) %in% subset_cell_IDs]
    }


    ## SCRAN ##
    marker_results <- scran::findMarkers(
        x = expr_data, groups = cell_metadata[[cluster_column]], ...
    )

    # data.table variables
    genes <- cluster <- feats <- NULL

    savelist <- lapply(names(marker_results), FUN = function(x) {
        dfr <- marker_results[[x]]
        DT <- data.table::as.data.table(dfr)
        DT[, feats := rownames(dfr)]
        DT[, cluster := x]
    })

    return(savelist)
}



#' @title findScranMarkers_one_vs_all
#' @name findScranMarkers_one_vs_all
#' @description Identify marker feats for all clusters in a one vs all manner
#' based on scran's implementation of findMarkers.
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values feat expression values to use
#' @param cluster_column clusters to use
#' @param subset_clusters subset of clusters to use
#' @param pval filter on minimal p-value
#' @param logFC filter on logFC
#' @param min_feats minimum feats to keep per cluster, overrides pval and logFC
#' @param min_genes deprecated, use min_feats
#' @param verbose be verbose
#' @param ...  additional parameters for the findMarkers function in scran
#' @returns data.table with marker feats
#' @seealso \code{\link{findScranMarkers}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findScranMarkers_one_vs_all(g, cluster_column = "leiden_clus")
#' @export
findScranMarkers_one_vs_all <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        cluster_column,
        subset_clusters = NULL,
        pval = 0.01,
        logFC = 0.5,
        min_feats = 10,
        min_genes = NULL,
        verbose = TRUE,
        ...) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
    }

    # verify if optional package is installed
    package_check(pkg_name = "scran", repository = "Bioc")

    # print message with information #
    if (verbose) {
        message("using 'Scran' to detect marker feats. If used in published
      research, please cite: Lun ATL, McCarthy DJ, Marioni JC (2016).
      'A step-by-step workflow for low-level analysis of single-cell RNA-seq
      data with Bioconductor.'
      F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2. ")
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

    # expression data
    values <- match.arg(
        expression_values,
        choices = unique(c(
            "normalized", "scaled", "custom",
            expression_values
        ))
    )

    # cluster column
    cell_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )
    if (!cluster_column %in% colnames(cell_metadata)) {
        stop("cluster column not found")
    }

    # restrict to a subset of clusters
    if (!is.null(subset_clusters)) {
        cell_metadata <- cell_metadata[get(cluster_column) %in% subset_clusters]
        subset_cell_IDs <- cell_metadata[["cell_ID"]]
        gobject <- subsetGiotto(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            cell_ids = subset_cell_IDs,
            verbose = FALSE
        )
        cell_metadata <- getCellMetadata(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "data.table",
            copy_obj = TRUE
        )
    }



    # sort uniq clusters
    uniq_clusters <- mixedsort(unique(cell_metadata[[cluster_column]]))


    # save list
    with_pbar({
        pb <- pbar(along = uniq_clusters)
        result_list <- lapply(
            seq_along(uniq_clusters),
            function(clus_i) {
                selected_clus <- uniq_clusters[clus_i]
                other_clus <- uniq_clusters[uniq_clusters != selected_clus]

                if (verbose == TRUE) {
                    cat("start with cluster ", selected_clus)
                }

                # one vs all markers
                markers <- findScranMarkers(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    expression_values = values,
                    cluster_column = cluster_column,
                    group_1 = selected_clus,
                    group_2 = other_clus,
                    verbose = FALSE
                )

                # identify list to continue with
                select_bool <- unlist(lapply(markers, FUN = function(x) {
                    unique(x$cluster) == selected_clus
                }))
                selected_table <- data.table::as.data.table(
                    markers[select_bool]
                )

                # remove summary column from scran output if present
                col_ind_keep <- !grepl("summary", colnames(selected_table))
                selected_table <- selected_table[, col_ind_keep, with = FALSE]

                # change logFC.xxx name to logFC
                data.table::setnames(
                    selected_table, colnames(selected_table)[4], "logFC"
                )
                data.table::setnames(
                    selected_table, colnames(selected_table)[5], "feats"
                )

                # filter selected table
                filtered_table <- selected_table[logFC > 0]
                filtered_table[, "ranking" := rank(-logFC)]

                # data.table variables
                p.value <- ranking <- NULL

                filtered_table <- filtered_table[
                    (p.value <= pval & logFC >= logFC) | (ranking <= min_feats)
                ]

                pb(message = c("cluster ", clus_i, "/", length(uniq_clusters)))
                return(filtered_table)
            }
        )
    })


    return(do.call("rbind", result_list))
}


#' @title findGiniMarkers
#' @name findGiniMarkers
#' @description Identify marker feats for selected clusters based on gini
#' detection and expression scores.
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values feat expression values to use
#' @param cluster_column clusters to use
#' @param subset_clusters selection of clusters to compare
#' @param group_1 group 1 cluster IDs from cluster_column for pairwise
#' comparison
#' @param group_1_name custom name for group_1 clusters
#' @param group_2 group 2 cluster IDs from cluster_column for pairwise
#' comparison
#' @param group_2_name custom name for group_2 clusters
#' @param min_expr_gini_score filter on minimum gini coefficient for expression
#' @param min_det_gini_score filter on minimum gini coefficient for detection
#' @param detection_threshold detection threshold for feat expression
#' @param rank_score rank scores for both detection and expression to include
#' @param min_feats minimum number of top feats to return
#' @param min_genes deprecated, use min_feats
#' @returns data.table with marker feats
#' @details
#' Detection of marker feats using the
#' [gini](https://en.wikipedia.org/wiki/Gini_coefficient)
#' coefficient is based on the following steps/principles per feat:
#' 1. calculate average expression per cluster
#' 2. calculate detection fraction per cluster
#' 3. calculate gini-coefficient for av. expression values over all clusters
#' 4. calculate gini-coefficient for detection fractions over all clusters
#' 5. convert gini-scores to rank scores
#' 6. for each feat create combined score = detection rank x expression rank x
#' expr gini-coefficient x detection gini-coefficient
#' 7. for each feat sort on expression and detection rank and combined score
#'
#' As a results "top gini" feats are feats that are very selectivily expressed
#' in a specific cluster,
#' however not always expressed in all cells of that cluster. In other words
#' highly specific, but
#' not necessarily sensitive at the single-cell level.
#'
#' To perform differential expression between custom selected groups of cells
#' you need to specify the cell_ID column to parameter \emph{cluster_column}
#' and provide the individual cell IDs to the parameters \emph{group_1} and
#' \emph{group_2}
#'
#' By default group names will be created by pasting the different id names
#' within each selected group.
#' When you have many different ids in a single group
#' it is recommend to provide names for both groups to \emph{group_1_name} and
#' \emph{group_2_name}
#' @md
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findGiniMarkers(g, cluster_column = "leiden_clus")
#' @export
findGiniMarkers <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        cluster_column,
        subset_clusters = NULL,
        group_1 = NULL,
        group_1_name = NULL,
        group_2 = NULL,
        group_2_name = NULL,
        min_expr_gini_score = 0.2,
        min_det_gini_score = 0.2,
        detection_threshold = 0,
        rank_score = 1,
        min_feats = 5,
        min_genes = NULL) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
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

    ## select expression values
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )


    # cluster column
    cell_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "cellMetaObj",
        copy_obj = TRUE
    )

    if (!cluster_column %in% colnames(cell_metadata[])) {
        stop("cluster column not found")
    }


    # subset clusters
    if (!is.null(subset_clusters)) {
        cell_metadata[] <- cell_metadata[][
            get(cluster_column) %in% subset_clusters
        ]
        subset_cell_IDs <- cell_metadata[][["cell_ID"]]
        gobject <- subsetGiotto(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            cell_ids = subset_cell_IDs
        )
    } else if (!is.null(group_1) & !is.null(group_2)) {
        cell_metadata[] <- cell_metadata[][
            get(cluster_column) %in% c(group_1, group_2)
        ]

        # create new pairwise group
        if (!is.null(group_1_name)) {
            if (!is.character(group_1_name)) {
                stop("group_1_name needs to be a character")
            }
            group_1_name <- group_1_name
        } else {
            group_1_name <- paste0(group_1, collapse = "_")
        }

        if (!is.null(group_2_name)) {
            if (!is.character(group_2_name)) {
                stop("group_2_name needs to be a character")
            }
            group_2_name <- group_2_name
        } else {
            group_2_name <- paste0(group_2, collapse = "_")
        }
        # data.table variables
        pairwise_select_comp <- NULL

        cell_metadata[][, pairwise_select_comp := ifelse(
            get(cluster_column) %in% group_1, group_1_name, group_2_name
        )]

        cluster_column <- "pairwise_select_comp"

        # expression data
        subset_cell_IDs <- cell_metadata[][["cell_ID"]]
        gobject <- subsetGiotto(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            cell_ids = subset_cell_IDs
        )

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- setGiotto(gobject, cell_metadata, verbose = FALSE)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    }


    # average expression per cluster
    aggr_sc_clusters <- create_average_DT(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        meta_data_name = cluster_column,
        expression_values = values
    )
    aggr_sc_clusters_DT <- data.table::as.data.table(aggr_sc_clusters)

    # data.table variables
    feats <- NULL

    aggr_sc_clusters_DT[, feats := rownames(aggr_sc_clusters)]
    aggr_sc_clusters_DT_melt <- data.table::melt.data.table(aggr_sc_clusters_DT,
        variable.name = "cluster",
        id.vars = "feats",
        value.name = "expression"
    )


    ## detection per cluster
    aggr_detection_sc_clusters <- create_average_detection_DT(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        meta_data_name = cluster_column,
        expression_values = values,
        detection_threshold = detection_threshold
    )
    aggr_detection_sc_clusters_DT <- data.table::as.data.table(
        aggr_detection_sc_clusters
    )
    aggr_detection_sc_clusters_DT[, feats := rownames(
        aggr_detection_sc_clusters
    )]
    aggr_detection_sc_clusters_DT_melt <- data.table::melt.data.table(
        aggr_detection_sc_clusters_DT,
        variable.name = "cluster",
        id.vars = "feats",
        value.name = "detection"
    )

    ## gini
    # data.table variables
    expression_gini <- detection_gini <- detection <- NULL

    aggr_sc_clusters_DT_melt[, expression_gini := mygini_fun(
        expression
    ), by = feats]
    aggr_detection_sc_clusters_DT_melt[, detection_gini := mygini_fun(
        detection
    ), by = feats]


    ## combine
    aggr_sc <- cbind(
        aggr_sc_clusters_DT_melt,
        aggr_detection_sc_clusters_DT_melt[
            , .(detection, detection_gini)
        ]
    )

    ## create combined rank

    # expression rank for each feat in all samples
    # rescale expression rank range between 1 and 0.1

    # data.table variables
    expression_rank <- cluster <- detection_rank <- NULL

    aggr_sc[, expression_rank := rank(-expression), by = feats]
    aggr_sc[, expression_rank := scales::rescale(
        expression_rank,
        to = c(1, 0.1)
    ), by = cluster]

    # detection rank for each feat in all samples
    # rescale detection rank range between 1 and 0.1
    aggr_sc[, detection_rank := rank(-detection), by = feats]
    aggr_sc[, detection_rank := scales::rescale(
        detection_rank,
        to = c(1, 0.1)
    ), by = cluster]

    # create combine score based on rescaled ranks and gini scores

    # data.table variables
    comb_score <- comb_rank <- NULL

    aggr_sc[, comb_score := (expression_gini * expression_rank) * (
        detection_gini * detection_rank)]
    setorder(aggr_sc, cluster, -comb_score)
    aggr_sc[, comb_rank := seq_len(.N), by = cluster]

    top_feats_scores <- aggr_sc[comb_rank <= min_feats | (
        expression_rank <= rank_score & detection_rank <= rank_score)]
    top_feats_scores_filtered <- top_feats_scores[comb_rank <= min_feats | (
        expression > min_expr_gini_score & detection > min_det_gini_score)]
    setorder(top_feats_scores_filtered, cluster, comb_rank)


    # remove 'cluster_' part if this is not part of the original cluster names
    original_uniq_cluster_names <- unique(cell_metadata[][[cluster_column]])
    if (sum(grepl("cluster_", original_uniq_cluster_names)) == 0) {
        top_feats_scores_filtered[, cluster := gsub(
            x = cluster, "cluster_", ""
        )]
    }

    return(top_feats_scores_filtered)
}




#' @title findGiniMarkers_one_vs_all
#' @name findGiniMarkers_one_vs_all
#' @description Identify marker feats for all clusters in a one vs all manner
#' based on gini detection and expression scores.
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values feat expression values to use
#' @param cluster_column clusters to use
#' @param subset_clusters selection of clusters to compare
#' @param min_expr_gini_score filter on minimum gini coefficient on expression
#' @param min_det_gini_score filter on minimum gini coefficient on detection
#' @param detection_threshold detection threshold for feat expression
#' @param rank_score rank scores for both detection and expression to include
#' @param min_feats minimum number of top feats to return
#' @param min_genes deprecated, use min_feats
#' @param verbose be verbose
#' @returns data.table with marker feats
#' @seealso \code{\link{findGiniMarkers}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findGiniMarkers_one_vs_all(g, cluster_column = "leiden_clus")
#' @export
findGiniMarkers_one_vs_all <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        cluster_column,
        subset_clusters = NULL,
        min_expr_gini_score = 0.5,
        min_det_gini_score = 0.5,
        detection_threshold = 0,
        rank_score = 1,
        min_feats = 4,
        min_genes = NULL,
        verbose = TRUE) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
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

    ## select expression values
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )


    # cluster column
    cell_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )

    if (!cluster_column %in% colnames(cell_metadata)) {
        stop("\n cluster column not found \n")
    }

    if (!is.null(subset_clusters)) {
        cell_metadata <- cell_metadata[get(cluster_column) %in% subset_clusters]
        subset_cell_IDs <- cell_metadata[["cell_ID"]]
        gobject <- subsetGiotto(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            cell_ids = subset_cell_IDs
        )
        cell_metadata <- getCellMetadata(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "data.table",
            copy_obj = TRUE
        )
    }


    # sort uniq clusters
    uniq_clusters <- mixedsort(unique(cell_metadata[[cluster_column]]))


    # GINI
    with_pbar({
        pb <- pbar(along = uniq_clusters)
        result_list <- lapply(
            seq_along(uniq_clusters),
            function(clus_i) {
                selected_clus <- uniq_clusters[clus_i]
                other_clus <- uniq_clusters[uniq_clusters != selected_clus]

                if (verbose == TRUE) {
                    cat("start with cluster ", selected_clus)
                }

                markers <- findGiniMarkers(
                    gobject = gobject,
                    feat_type = feat_type,
                    spat_unit = spat_unit,
                    expression_values = values,
                    cluster_column = cluster_column,
                    group_1 = selected_clus,
                    group_2 = other_clus,
                    min_expr_gini_score = min_expr_gini_score,
                    min_det_gini_score = min_det_gini_score,
                    detection_threshold = detection_threshold,
                    rank_score = rank_score,
                    min_feats = min_feats
                )

                # filter steps

                # data.table variables
                cluster <- NULL

                filtered_table <- markers[cluster == selected_clus]

                pb(message = c("cluster ", clus_i, "/", length(uniq_clusters)))
                return(filtered_table)
            }
        )
    })

    return(do.call("rbind", result_list))
}



#' @title findMastMarkers
#' @name findMastMarkers
#' @description Identify marker feats for selected clusters based on the
#' MAST package.
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values feat expression values to use
#' @param cluster_column clusters to use
#' @param group_1 group 1 cluster IDs from cluster_column for pairwise
#' comparison
#' @param group_1_name custom name for group_1 clusters
#' @param group_2 group 2 cluster IDs from cluster_column for pairwise
#' comparison
#' @param group_2_name custom name for group_2 clusters
#' @param adjust_columns column in pDataDT to adjust for (e.g. detection rate)
#' @param verbose be verbose
#' @param ... additional parameters for the zlm function in MAST
#' @return data.table with marker feats
#' @details This is a minimal convenience wrapper around the
#' \code{\link[MAST]{zlm}}
#' from the MAST package to detect differentially expressed feats. Caution:
#' with large datasets
#' MAST might take a long time to run and finish
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findMastMarkers(
#'     gobject = g, cluster_column = "leiden_clus", group_1 = 1,
#'     group_2 = 2
#' )
#' @export
findMastMarkers <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        cluster_column,
        group_1 = NULL,
        group_1_name = NULL,
        group_2 = NULL,
        group_2_name = NULL,
        adjust_columns = NULL,
        verbose = FALSE,
        ...) {
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

    # verify if optional package is installed
    package_check(pkg_name = "MAST", repository = "Bioc")

    # print message with information #
    if (verbose) {
        message("using 'MAST' to detect marker feats. If used in published
      research, please cite: McDavid A, Finak G, Yajima M (2020).
      MAST: Model-based Analysis of Single Cell Transcriptomics.
      R package version 1.14.0, https://github.com/RGLab/MAST/.")
    }

    ## select expression values to use
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )

    ## cluster column
    cell_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "cellMetaObj",
        copy_obj = TRUE
    )
    if (!cluster_column %in% colnames(cell_metadata[])) {
        stop("cluster column not found")
    }

    ## select group ids
    if (is.null(group_1) | is.null(group_2)) {
        stop("specificy group ids for both group_1 and group_2")
    }

    ## subset data based on group_1 and group_2
    cell_metadata[] <- cell_metadata[][
        get(cluster_column) %in% c(group_1, group_2)
    ]
    if (nrow(cell_metadata[]) == 0) {
        stop("there are no cells for group_1 or group_2, check cluster column")
    }

    ## create new pairwise group
    if (is.null(group_1_name)) group_1_name <- paste0(group_1, collapse = "_")
    if (is.null(group_2_name)) group_2_name <- paste0(group_2, collapse = "_")

    # data.table variables
    pairwise_select_comp <- NULL

    cell_metadata[][, pairwise_select_comp := ifelse(
        get(cluster_column) %in% group_1, group_1_name, group_2_name
    )]

    if (nrow(cell_metadata[][pairwise_select_comp == group_1_name]) == 0) {
        stop("there are no cells for group_1, check cluster column")
    }

    if (nrow(cell_metadata[][pairwise_select_comp == group_2_name]) == 0) {
        stop("there are no cells for group_2, check cluster column")
    }

    cluster_column <- "pairwise_select_comp"

    # expression data
    subset_cell_IDs <- cell_metadata[][["cell_ID"]]
    gobject <- subsetGiotto(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        cell_ids = subset_cell_IDs
    )

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject <- setGiotto(gobject, cell_metadata, verbose = FALSE)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


    ## START MAST ##
    ## create mast object ##
    # expression data
    values <- match.arg(
        expression_values,
        choices = unique(c(
            "normalized", "scaled", "custom",
            expression_values
        ))
    )
    expr_data <- getExpression(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        values = values,
        output = "matrix"
    )
    # column & row data
    column_data <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )
    setnames(column_data, "cell_ID", "wellKey")
    row_data <- getFeatureMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )
    setnames(row_data, "feat_ID", "primerid")
    # mast object
    mast_data <- MAST::FromMatrix(
        exprsArray = as.matrix(expr_data)[, column_data[["wellKey"]]],
        cData = column_data,
        fData = row_data,
        check_sanity = FALSE
    )

    ## set conditions and relevel
    cond <- factor(SingleCellExperiment::colData(mast_data)[[cluster_column]])
    cond <- stats::relevel(cond, group_2_name)
    mast_data@colData[[cluster_column]] <- cond

    ## create formula and run MAST feat regressions
    if (!is.null(adjust_columns)) {
        myformula <- stats::as.formula(paste0(
            "~ 1 + ", cluster_column, " + ",
            paste(adjust_columns, collapse = " + ")
        ))
    } else {
        myformula <- stats::as.formula(paste0("~ 1 + ", cluster_column))
    }
    zlmCond <- MAST::zlm(formula = myformula, sca = mast_data, ...)

    ## run LRT and return data.table with results

    # data.table variables
    contrast <- component <- primerid <- `Pr(>Chisq)` <- coef <-
        ci.hi <- ci.lo <- fdr <- NULL

    sample <- paste0(cluster_column, group_1_name)
    summaryCond <- MAST::summary(zlmCond, doLRT = sample)
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(
        summaryDt[
            contrast == sample & component == "H",
            .(primerid, `Pr(>Chisq)`)
        ], # hurdle P values
        summaryDt[
            contrast == sample & component == "logFC",
            .(primerid, coef, ci.hi, ci.lo)
        ],
        by = "primerid"
    ) # logFC coefficients
    fcHurdle[, fdr := stats::p.adjust(`Pr(>Chisq)`, "fdr")]
    data.table::setorder(fcHurdle, fdr)

    # data.table variables
    cluster <- NULL

    fcHurdle[, cluster := paste0(group_1_name, "_vs_", group_2_name)]
    data.table::setnames(fcHurdle, old = "primerid", new = "feats")

    return(fcHurdle)
}




#' @title findMastMarkers_one_vs_all
#' @name findMastMarkers_one_vs_all
#' @description Identify marker feats for all clusters in a one vs all manner
#' based on the MAST package.
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values feat expression values to use
#' @param cluster_column clusters to use
#' @param subset_clusters selection of clusters to compare
#' @param adjust_columns column in pDataDT to adjust for (e.g. detection rate)
#' @param pval filter on minimal p-value
#' @param logFC filter on logFC
#' @param min_feats minimum feats to keep per cluster, overrides pval and logFC
#' @param min_genes deprecated, use min_feats
#' @param verbose be verbose
#' @param ... additional parameters for the zlm function in MAST
#' @returns data.table with marker feats
#' @seealso \code{\link{findMastMarkers}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findMastMarkers_one_vs_all(gobject = g, cluster_column = "leiden_clus")
#' @export
findMastMarkers_one_vs_all <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        cluster_column,
        subset_clusters = NULL,
        adjust_columns = NULL,
        pval = 0.001,
        logFC = 1,
        min_feats = 10,
        min_genes = NULL,
        verbose = TRUE,
        ...) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
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

    # verify if optional package is installed
    package_check(pkg_name = "MAST", repository = "Bioc")

    # print message with information #
    if (verbose) {
        message("using 'MAST' to detect marker feats. If used in published
        research, please cite: McDavid A, Finak G, Yajima M (2020).
        MAST: Model-based Analysis of Single Cell Transcriptomics.
        R package version 1.14.0, https://github.com/RGLab/MAST/.")
    }


    ## cluster column
    cell_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )
    if (!cluster_column %in% colnames(cell_metadata)) {
        stop("cluster column not found")
    }

    # restrict to a subset of clusters
    if (!is.null(subset_clusters)) {
        cell_metadata <- cell_metadata[get(cluster_column) %in% subset_clusters]
        subset_cell_IDs <- cell_metadata[["cell_ID"]]
        gobject <- subsetGiotto(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            cell_ids = subset_cell_IDs,
            verbose = FALSE
        )
        cell_metadata <- getCellMetadata(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            output = "data.table",
            copy_obj = TRUE
        )
    }

    ## sort uniq clusters
    uniq_clusters <- mixedsort(unique(cell_metadata[[cluster_column]]))

    # save list
    result_list <- list()

    for (clus_i in seq_along(uniq_clusters)) {
        selected_clus <- uniq_clusters[clus_i]
        other_clus <- uniq_clusters[uniq_clusters != selected_clus]

        if (verbose == TRUE) {
            cat("start with cluster ", selected_clus)
        }

        temp_mast_markers <- findMastMarkers(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            expression_values = expression_values,
            cluster_column = cluster_column,
            adjust_columns = adjust_columns,
            group_1 = selected_clus,
            group_1_name = selected_clus,
            group_2 = other_clus,
            group_2_name = "others",
            verbose = FALSE
        )

        result_list[[clus_i]] <- temp_mast_markers
    }

    # filter or retain only selected marker feats
    result_dt <- do.call("rbind", result_list)

    # data.table variables
    ranking <- fdr <- coef <- NULL

    result_dt[, ranking := seq_len(.N), by = "cluster"]
    filtered_result_dt <- result_dt[
        ranking <= min_feats | (fdr < pval & coef > logFC)
    ]

    return(filtered_result_dt)
}






#' @title findMarkers
#' @name findMarkers
#' @description Identify marker feats for selected clusters.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param expression_values feat expression values to use
#' @param cluster_column clusters to use
#' @param method method to use to detect differentially expressed feats
#' @param subset_clusters selection of clusters to compare
#' @param group_1 group 1 cluster IDs from cluster_column for pairwise
#' comparison
#' @param group_2 group 2 cluster IDs from cluster_column for pairwise
#' comparison
#' @param min_expr_gini_score gini: filter on minimum gini coefficient for
#' expression
#' @param min_det_gini_score gini: filter minimum gini coefficient for detection
#' @param detection_threshold gini: detection threshold for feat expression
#' @param rank_score gini: rank scores to include
#' @param min_feats minimum number of top feats to return (for gini)
#' @param min_genes deprecated, use min_feats
#' @param group_1_name mast: custom name for group_1 clusters
#' @param group_2_name mast: custom name for group_2 clusters
#' @param adjust_columns mast: column in pDataDT to adjust for
#' (e.g. detection rate)
#' @param ... additional parameters for the findMarkers function in scran or
#' zlm function in MAST
#' @returns data.table with marker feats
#' @details Wrapper for all individual functions to detect marker feats for
#' clusters.
#' @seealso \code{\link{findScranMarkers}}, \code{\link{findGiniMarkers}} and
#' \code{\link{findMastMarkers}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findMarkers(g, cluster_column = "leiden_clus")
#' @export
findMarkers <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        cluster_column = NULL,
        method = c("scran", "gini", "mast"),
        subset_clusters = NULL,
        group_1 = NULL,
        group_2 = NULL,
        min_expr_gini_score = 0.5,
        min_det_gini_score = 0.5,
        detection_threshold = 0,
        rank_score = 1,
        min_feats = 4,
        min_genes = NULL,
        group_1_name = NULL,
        group_2_name = NULL,
        adjust_columns = NULL,
        ...) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
    }

    # input
    if (is.null(cluster_column)) {
        stop("A valid cluster column needs to be given to cluster_column,
            see pDataDT()")
    }

    # select method
    method <- match.arg(method, choices = c("scran", "gini", "mast"))

    if (method == "scran") {
        markers_result <- findScranMarkers(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            expression_values = expression_values,
            cluster_column = cluster_column,
            subset_clusters = subset_clusters,
            group_1 = group_1,
            group_2 = group_2,
            group_1_name = group_1_name,
            group_2_name = group_2_name,
            ...
        )
    } else if (method == "gini") {
        markers_result <- findGiniMarkers(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            expression_values = expression_values,
            cluster_column = cluster_column,
            subset_clusters = subset_clusters,
            group_1 = group_1,
            group_2 = group_2,
            group_1_name = group_1_name,
            group_2_name = group_2_name,
            min_expr_gini_score = min_expr_gini_score,
            min_det_gini_score = min_det_gini_score,
            detection_threshold = detection_threshold,
            rank_score = rank_score,
            min_feats = min_feats
        )
    } else if (method == "mast") {
        markers_result <- findMastMarkers(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            expression_values = expression_values,
            cluster_column = cluster_column,
            group_1 = group_1,
            group_1_name = group_1_name,
            group_2 = group_2,
            group_2_name = group_2_name,
            adjust_columns = adjust_columns,
            ...
        )
    }

    return(markers_result)
}


#' @title findMarkers_one_vs_all
#' @name findMarkers_one_vs_all
#' @description Identify marker feats for all clusters in a one vs all manner.
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values feat expression values to use
#' @param cluster_column clusters to use
#' @param method method to use to detect differentially expressed feats
#' @param subset_clusters selection of clusters to compare
#' @param pval scran & mast: filter on minimal p-value
#' @param logFC scan & mast: filter on logFC
#' @param min_feats minimum feats to keep per cluster, overrides pval and logFC
#' @param min_genes deprecated, use min_feats
#' @param min_expr_gini_score gini: filter on minimum gini coefficient for
#' expression
#' @param min_det_gini_score gini: filter minimum gini coefficient for detection
#' @param detection_threshold gini: detection threshold for feat expression
#' @param rank_score gini: rank scores to include
#' @param adjust_columns mast: column in pDataDT to adjust for
#' (e.g. detection rate)
#' @param verbose be verbose
#' @param ... additional parameters for the findMarkers function in scran or
#' zlm function in MAST
#' @returns data.table with marker feats
#' @details Wrapper for all one vs all functions to detect marker feats for
#' clusters.
#' @seealso \code{\link{findScranMarkers_one_vs_all}},
#' \code{\link{findGiniMarkers_one_vs_all}} and
#' \code{\link{findMastMarkers_one_vs_all}}
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' findMarkers_one_vs_all(g, cluster_column = "leiden_clus")
#' @export
findMarkers_one_vs_all <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        cluster_column,
        subset_clusters = NULL,
        method = c("scran", "gini", "mast"),
        # scran & mast
        pval = 0.01,
        logFC = 0.5,
        min_feats = 10,
        min_genes = NULL,
        # gini
        min_expr_gini_score = 0.5,
        min_det_gini_score = 0.5,
        detection_threshold = 0,
        rank_score = 1,
        # mast specific
        adjust_columns = NULL,
        verbose = TRUE,
        ...) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
    }

    # select method
    method <- match.arg(method, choices = c("scran", "gini", "mast"))

    if (method == "scran") {
        markers_result <- findScranMarkers_one_vs_all(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            expression_values = expression_values,
            cluster_column = cluster_column,
            subset_clusters = subset_clusters,
            pval = pval,
            logFC = logFC,
            min_feats = min_feats,
            verbose = verbose,
            ...
        )
    } else if (method == "gini") {
        markers_result <- findGiniMarkers_one_vs_all(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            expression_values = expression_values,
            cluster_column = cluster_column,
            subset_clusters = subset_clusters,
            min_expr_gini_score = min_expr_gini_score,
            min_det_gini_score = min_det_gini_score,
            detection_threshold = detection_threshold,
            min_feats = min_feats,
            verbose = verbose
        )
    } else if (method == "mast") {
        markers_result <- findMastMarkers_one_vs_all(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            expression_values = expression_values,
            cluster_column = cluster_column,
            subset_clusters = subset_clusters,
            adjust_columns = adjust_columns,
            pval = pval,
            logFC = logFC,
            min_feats = min_feats,
            verbose = verbose,
            ...
        )
    }

    return(markers_result)
}
