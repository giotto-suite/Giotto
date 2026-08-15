# markersParam ####

#' @name markers_scran
#' @title Pairwise Marker Detection (scran)
#' @description
#' Detect marker features by comparing each group against the others, as
#' implemented by \code{\link[scran]{findMarkers}}.
#'
#' Every ordered pair of groups is tested, and the resulting p-values are
#' combined into a single ranked table per group. With the default
#' `test_type = "t"` the test is a Welch \eqn{t}-test on the per-group
#' moments:
#'
#' \deqn{\LARGE
#' t = \frac{\bar{x}_{i,a} - \bar{x}_{i,b}}
#'          {\sqrt{s^2_{i,a}/n_a + s^2_{i,b}/n_b}}
#' }
#' Where:
#'
#' * (\eqn{\bar{x}_{i,a}}) is the mean of feature \eqn{i} over the cells of
#' group \eqn{a}
#' * (\eqn{s^2_{i,a}}) is the variance of feature \eqn{i} over the same cells
#' * (\eqn{n_a}) is the number of cells in group \eqn{a}
#'
#' Because the statistic depends on the values only through
#' \eqn{n}, \eqn{\bar{x}} and \eqn{s^2}, the expression matrix is visited
#' once per analysis rather than once per comparison. That is what lets a
#' streaming backend implement the same test without materializing the matrix.
#'
#' Which expression values are tested is the caller's choice — the methods use
#' whatever matrix they are given.
#' @section params:
#'
#' \tabular{ll}{
#'   `test_type` \tab character (default = "t"). Pairwise test: `"t"` Welch
#'   \eqn{t}-test, `"wilcox"` rank-sum, `"binom"` binomial on detection
#'   rates. Backends may support only a subset. \cr
#'   `pval_type` \tab character (default = "any"). How the pairwise p-values
#'   are combined per group: `"any"`, `"some"`, `"all"`. \cr
#'   `comparison` \tab character (default = "pairwise"). `"pairwise"` tests
#'   every ordered pair; `"one_vs_rest"` tests each group against the pooled
#'   remainder. \cr
#'   `direction` \tab character (default = "any"). `"any"`, `"up"`, `"down"`.
#'   \cr
#'   `lfc` \tab numeric (default = 0). Log-fold-change threshold to test
#'   against. \cr
#'   `std_lfc` \tab logical (default = FALSE). Report the effect size as a
#'   standardized log-fold-change (Cohen's d). \cr
#'   `min_prop` \tab numeric or NULL. Minimum proportion of comparisons a
#'   feature must be significant in, for `pval_type = "some"`. \cr
#'   `log_p` \tab logical (default = FALSE). Report p-values on the log
#'   scale. \cr
#'   `full_stats` \tab logical (default = FALSE). Retain the per-comparison
#'   statistics as nested columns. \cr
#'   `sorted` \tab logical (default = TRUE). Sort each group's table by
#'   significance.
#' }
#'
#' `pval_type` and `min_prop` describe how scran combines the pairwise
#' comparisons and have no meaning independent of it; see
#' \code{\link[scran]{findMarkers}}.
#' @md
#' @family marker detection parameters
#' @seealso [analyze_param], [markersParam()], [findScranMarkers()]
#' @returns marker detection results
NULL


#' @rdname analyze_param
#' @exportClass markersParam
setClass("markersParam", contains = c("VIRTUAL", "analyzeParam"))

#' @rdname analyze_param
#' @exportClass scranMarkersParam
setClass("scranMarkersParam", contains = "markersParam")


# param factory ####

#' @rdname analyze_param
#' @export
markersParam <- function(method = "scran", ...) {
    method <- match.arg(tolower(method), c("scran"))
    switch(method,
        "scran" = .markers_param_scran(...)
    )
}

#' @keywords internal
#' @noRd
.markers_param_scran <- function(...) {
    p <- new("scranMarkersParam", param = list(...))
    p$test_type <- p$test_type %null% "t"
    p$pval_type <- p$pval_type %null% "any"
    p$comparison <- p$comparison %null% "pairwise"
    p$direction <- p$direction %null% "any"
    p$lfc <- as.numeric(p$lfc %null% 0)
    p$std_lfc <- isTRUE(p$std_lfc)
    p$log_p <- isTRUE(p$log_p)
    p$full_stats <- isTRUE(p$full_stats)
    p$sorted <- p$sorted %null% TRUE
    p
}


# analyzeData(<matrix>, scranMarkersParam) ####

# In-memory marker detection is a pass-through to scran, deliberately.
#
# scran's `findMarkers(test.type = "t")` already makes ONE optimal C++ pass
# over the matrix for its per-group moments, so there is nothing to gain by
# reimplementing it here -- and a great deal to lose, since a second
# implementation would have to be kept in step with scran's forever. Backends
# that cannot hand scran a matrix at all (streaming stores) carry their own
# equivalent; this one must not.
#
#' @rdname markers_scran
#' @param x expression values. A `matrix`, a `Matrix`, or a `DelayedMatrix`
#'   (which is what `expression_values = "scaled"` holds).
#' @param param a [scranMarkersParam-class].
#' @param groups vector of group assignments, one per cell, in column order.
#' @export
setMethod("analyzeData",
    signature(x = "matrix", param = "scranMarkersParam"),
    function(x, param, ..., groups = NULL) {
        .markers_scran(x, param, groups = groups)
    }
)

#' @rdname markers_scran
#' @export
setMethod("analyzeData",
    signature(x = "Matrix", param = "scranMarkersParam"),
    function(x, param, ..., groups = NULL) {
        .markers_scran(x, param, groups = groups)
    }
)

# Covers `ScaledMatrix`, which is what `expression_values = "scaled"` holds.
#' @rdname markers_scran
#' @export
setMethod("analyzeData",
    signature(x = "DelayedMatrix", param = "scranMarkersParam"),
    function(x, param, ..., groups = NULL) {
        .markers_scran(x, param, groups = groups)
    }
)


#' @keywords internal
#' @noRd
.markers_scran <- function(x, param, groups) {
    package_check(pkg_name = "scran", repository = "Bioc")
    if (is.null(groups)) {
        stop("[analyzeData(scranMarkersParam)] `groups` is required: one group ",
            "assignment per cell.", call. = FALSE)
    }
    if (identical(param$comparison %null% "pairwise", "one_vs_rest")) {
        return(.markers_one_vs_rest_scran(x, groups, param))
    }
    do.call(scran::findMarkers,
        c(list(x = x, groups = groups), .markers_scran_args(param)))
}


# Translate the param's snake_case fields to scran's dotted argument names.
# Anything else on `@param` is passed through untouched, which is what lets
# `findScranMarkers(...)` keep forwarding arbitrary scran arguments.
#' @keywords internal
#' @noRd
.markers_scran_args <- function(param) {
    p <- as.list(param@param)
    rename <- c(
        test_type = "test.type", pval_type = "pval.type",
        std_lfc = "std.lfc", min_prop = "min.prop",
        log_p = "log.p", full_stats = "full.stats"
    )
    # Not scran arguments: `comparison` selects which sweep this method runs.
    p[["comparison"]] <- NULL
    for (from in names(rename)) {
        if (from %in% names(p)) {
            p[[rename[[from]]]] <- p[[from]]
            p[[from]] <- NULL
        }
    }
    p[!vapply(p, is.null, logical(1L))]
}


# One scran call per group, each against the pooled remainder.
#
# This is G accumulator passes and stays that way: pooling the moments would
# be one pass, but scran exposes no way to inject precomputed moments, so
# sharing the pass in memory would mean transcribing its statistics. Level
# names match the streaming backend so both produce the same table shape.
#' @keywords internal
#' @noRd
.markers_one_vs_rest_scran <- function(x, groups, param) {
    lvls <- levels(droplevels(
        if (is.factor(groups)) groups else factor(groups)
    ))
    args <- .markers_scran_args(param)
    out <- lapply(lvls, function(k) {
        rest <- setdiff(lvls, k)
        pooled <- ifelse(groups == k, k, paste0(rest, collapse = "_"))
        res <- do.call(scran::findMarkers,
            c(list(x = x, groups = pooled), args))
        res[[k]]
    })
    names(out) <- lvls
    S4Vectors::SimpleList(out)
}


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


    ## MARKERS ##
    # Dispatched on the expression object, not called directly: `getExpression`
    # returns whatever the slot holds, so a disk-backed store arrives here
    # intact and routes to its own streaming method. In memory this is a
    # pass-through to `scran::findMarkers`.
    marker_results <- analyzeData(
        x = expr_data,
        param = markersParam(method = "scran", ...),
        groups = cell_metadata[[cluster_column]]
    )

    lapply(names(marker_results), function(x) {
        .markers_result_dt(marker_results[[x]], cluster = x)
    })
}


# The one place a marker DataFrame becomes a data.table, shared by
# `findScranMarkers()` and `findScranMarkers_one_vs_all()`. Both backends
# return scran's shape, so neither needs to know which produced it.
#' @keywords internal
#' @noRd
.markers_result_dt <- function(dfr, cluster) {
    # data.table variables
    feats <- NULL

    DT <- data.table::as.data.table(dfr)
    DT[, feats := rownames(dfr)]
    DT[, "cluster" := cluster]
    DT[]
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


    # One dispatched call for the whole sweep, rather than one per cluster.
    #
    # The comparison is unchanged -- each cluster against the pooled remainder
    # -- but the backend now decides how to get there. A streaming store takes
    # one statistic pass and pools the rest group arithmetically; the in-memory
    # method still runs one scran call per cluster, because scran offers no way
    # to inject precomputed moments and sharing the pass there would mean
    # reimplementing its statistics.
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )
    marker_results <- analyzeData(
        x = expr_data,
        param = markersParam(method = "scran", comparison = "one_vs_rest"),
        groups = cell_metadata[[cluster_column]]
    )

    # save list
    with_pbar({
        pb <- pbar(along = uniq_clusters)
        result_list <- lapply(
            seq_along(uniq_clusters),
            function(clus_i) {
                selected_clus <- uniq_clusters[clus_i]

                if (verbose == TRUE) {
                    cat("start with cluster ", selected_clus)
                }

                selected_table <- .markers_result_dt(
                    marker_results[[as.character(selected_clus)]],
                    cluster = selected_clus
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
                selected_table[, "ranking" := rank(-logFC)]

                # data.table variables
                p.value <- ranking <- NULL

                selected_table <- selected_table[
                    (p.value <= pval & logFC >= logFC) | (ranking <= min_feats)
                ]

                pb(message = c("cluster ", clus_i, "/", length(uniq_clusters)))
                return(selected_table)
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
