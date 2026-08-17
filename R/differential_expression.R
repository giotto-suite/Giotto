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


#' @name markers_gini
#' @title Specificity Marker Detection (gini)
#' @description
#' Detect marker features by how unevenly a feature is distributed across
#' groups, using the Gini coefficient.
#'
#' Two coefficients are taken per feature, over the per-group mean expression
#' and over the per-group detection fraction:
#'
#' \deqn{\LARGE
#' G_i = \frac{\sum_{a}\sum_{b} |v_{i,a} - v_{i,b}|}
#'            {2 G^2 \bar{v}_i}
#' }
#' Where:
#'
#' * (\eqn{v_{i,a}}) is the statistic for feature \eqn{i} in group \eqn{a} —
#' mean expression for `expression_gini`, detection fraction for
#' `detection_gini`
#' * (\eqn{G}) is the number of groups
#'
#' The statistic depends on the values only through the per-(feature, group)
#' mean and detection fraction, so the expression matrix is visited once
#' regardless of how many groups there are. That is what lets a streaming
#' backend supply gini markers with no code of its own — the pass is
#' [featStatsParam-class], and this method consumes its output.
#'
#' Gini is **scale-free**: a 0.001 vs 0.0001 difference between groups scores
#' identically to 100 vs 10. It therefore carries no magnitude term and will
#' rank near-noise features as perfectly specific, which is what
#' `min_expression` and `min_detection` exist to prevent. Its ceiling is
#' \eqn{(G-1)/G}, so a fixed `min_expression_gini` is not comparable across
#' runs with different group counts unless `min_length` is set.
#'
#' Which expression values are scored is the caller's choice — the method uses
#' whatever matrix it is given.
#' @section params:
#'
#' \tabular{ll}{
#'   `comparison` \tab character (default = "pairwise"). `"pairwise"` scores
#'   every group against the others at once; `"one_vs_rest"` scores each group
#'   against the pooled remainder, one table per group. \cr
#'   `min_expression` \tab numeric (default = 0.2). Minimum per-group mean
#'   expression, gating the `expression` column. \cr
#'   `min_detection` \tab numeric (default = 0.2). Minimum fraction of a
#'   group's cells above `detection_threshold`, gating `detection`. \cr
#'   `min_expression_gini` \tab numeric (default = -Inf). Minimum coefficient,
#'   gating `expression_gini`. \cr
#'   `min_detection_gini` \tab numeric (default = -Inf). Minimum coefficient,
#'   gating `detection_gini`. \cr
#'   `detection_threshold` \tab numeric (default = 0). Value above which a cell
#'   counts as expressing a feature. Not a filter on returned rows. \cr
#'   `min_length` \tab integer (default = 0). Pad each per-group vector to this
#'   length with copies of its minimum before taking the coefficient, removing
#'   the dependence on group count so scores compare across runs. `0` never
#'   pads. \cr
#'   `rank_score` \tab numeric (default = Inf). Keep a feature when its group is
#'   within this rank for both expression and detection. \cr
#'   `min_feats` \tab integer (default = 5). Keep this many top features per
#'   group regardless of the gates, so a group is never empty.
#' }
#'
#' The gates are OR'd with `min_feats`, so tightening them shrinks the result
#' toward `min_feats` per group and never below it.
#'
#' Defaults here are [findGiniMarkers()]'s. [findGiniMarkers_one_vs_all()]
#' passes `min_expression = 0.5`, `min_detection = 0.5` and `min_feats = 4`
#' explicitly; `comparison = "one_vs_rest"` does **not** switch them for you.
#' @md
#' @family marker detection parameters
#' @seealso [analyze_param], [markersParam()], [findGiniMarkers()],
#'   [featStatsParam-class]
#' @returns marker detection results
NULL


#' @rdname analyze_param
#' @exportClass markersParam
setClass("markersParam", contains = c("VIRTUAL", "analyzeParam"))

#' @rdname analyze_param
#' @exportClass scranMarkersParam
setClass("scranMarkersParam", contains = "markersParam")

#' @rdname analyze_param
#' @exportClass giniMarkersParam
setClass("giniMarkersParam", contains = "markersParam")


# param factory ####

#' @rdname analyze_param
#' @export
markersParam <- function(method = "scran", ...) {
    method <- match.arg(tolower(method), c("scran", "gini"))
    switch(method,
        "scran" = .markers_param_scran(...),
        "gini" = .markers_param_gini(...)
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

#' @keywords internal
#' @noRd
.markers_param_gini <- function(...) {
    p <- new("giniMarkersParam", param = list(...))
    p$comparison <- p$comparison %null% "pairwise"
    p$min_expression <- as.numeric(p$min_expression %null% 0.2)
    p$min_detection <- as.numeric(p$min_detection %null% 0.2)
    p$min_expression_gini <- as.numeric(p$min_expression_gini %null% -Inf)
    p$min_detection_gini <- as.numeric(p$min_detection_gini %null% -Inf)
    p$detection_threshold <- as.numeric(p$detection_threshold %null% 0)
    p$min_length <- as.numeric(p$min_length %null% 0)
    p$rank_score <- as.numeric(p$rank_score %null% Inf)
    p$min_feats <- as.numeric(p$min_feats %null% 5)
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
#' @param groups vector of group assignments, one per cell. Taken in column
#'   order of `x`, unless the vector is named by cell ID, in which case it is
#'   matched to `colnames(x)`. Naming is the safer form: the expression matrix
#'   and the cell metadata a caller builds `groups` from are fetched
#'   independently and need not share a cell order.
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
    # positional against the columns of `x` unless named by cell ID
    groups <- .align_groups(x, groups)

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

    # `std.lfc` is t-test only -- `pairwiseWilcox()` and `pairwiseBinom()` do
    # not take it, and an effect size expressed in pooled standard deviations
    # has no counterpart for an AUC or a detection rate. The param always
    # carries the field (defaulted to FALSE), so forwarding it unconditionally
    # breaks every non-t test with "unused argument (std.lfc = FALSE)".
    if (!identical(p$test_type %null% "t", "t")) {
        if (isTRUE(p$std_lfc)) {
            stop("[markersParam] `std_lfc = TRUE` applies only to ",
                "`test_type = \"t\"`; ", p$test_type, " reports no ",
                "standardized effect size.", call. = FALSE)
        }
        p[["std_lfc"]] <- NULL
    }

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
        # named so the verb matches on cell ID rather than position
        groups = stats::setNames(
            cell_metadata[[cluster_column]], cell_metadata[["cell_ID"]]
        )
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
        groups = stats::setNames(
            cell_metadata[[cluster_column]], cell_metadata[["cell_ID"]]
        )
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

                # `logFC` names both this function's argument and a column of
                # `selected_table`, and inside `i` the column wins -- so
                # `logFC >= logFC` compared the column against itself and was
                # always TRUE. Bind the threshold to a name that cannot
                # collide.
                logFC_thresh <- logFC

                selected_table <- selected_table[
                    (p.value <= pval & logFC >= logFC_thresh) |
                        (ranking <= min_feats)
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
#' @param min_expression minimum per-cluster mean expression, gating the
#' `expression` column of the result
#' @param min_detection minimum fraction of a cluster's cells with expression
#' above `detection_threshold`, gating the `detection` column of the result
#' @param min_expression_gini minimum gini coefficient of expression, gating
#' the `expression_gini` column of the result. `-Inf` (default) disables it.
#' @param min_detection_gini minimum gini coefficient of detection, gating the
#' `detection_gini` column of the result. `-Inf` (default) disables it.
#' @param detection_threshold expression value above which a cell counts as
#' expressing a feature, used when computing `detection`. Not a filter on the
#' returned rows -- see `min_detection` for that.
#' @param min_length pad the per-cluster vector to this length before taking
#' the gini coefficient, using copies of its minimum. Removes the dependence
#' of the coefficient on how many clusters were compared, so gini scores and
#' the `min_expression_gini` / `min_detection_gini` thresholds become
#' comparable across runs. `0` (the default) never pads.
#' @param rank_score keep a feature when its cluster is within this rank for
#' both `expression` and `detection`, where rank 1 is the cluster in which the
#' feature is highest. `Inf` (default) disables it. Combined with `min_feats`
#' by `or`, like the other gates.
#' @param min_feats minimum number of top feats to return
#' @param min_genes deprecated, use min_feats
#' @param min_expr_gini_score `r lifecycle::badge("deprecated")` use
#' `min_expression`. Despite its name it never gated a gini coefficient.
#' @param min_det_gini_score `r lifecycle::badge("deprecated")` use
#' `min_detection`. Despite its name it never gated a gini coefficient.
#' @returns data.table with marker feats
#' @details
#' Detection of marker feats using the
#' [gini](https://en.wikipedia.org/wiki/Gini_coefficient)
#' coefficient is based on the following steps/principles, per feat:
#' 1. `expression`: the mean expression over the cells of each cluster.
#' 2. `detection`: the fraction of each cluster's cells whose expression
#'    exceeds `detection_threshold`.
#' 3. `expression_gini`: the gini coefficient of the `expression` values.
#' 4. `detection_gini`: the gini coefficient of the `detection` values.
#' 5. `expression_rank` / `detection_rank`: the clusters ranked by `expression`
#'    and by `detection`, rank 1 being the cluster where the feature is
#'    highest. These are ranks of the **values**, not of the gini
#'    coefficients, and are what `rank_score` filters on. Ties share a rank.
#' 6. `comb_score`: the same two orderings are rescaled within each cluster to
#'    \[1, 0.1\] to act as weights, and `comb_score` = `expression_gini` x
#'    expression weight x `detection_gini` x detection weight. The weights are
#'    internal and are not the returned rank columns: they resolve ties by
#'    averaging where step 5 resolves them by minimum, so on a tied feature
#'    `comb_score` cannot be reproduced from the returned ranks alone.
#' 7. within each cluster, sort by `comb_score` (descending) and number the
#'    result as `comb_rank`.
#'
#' Steps 1 and 2 reduce the matrix to two features x clusters tables; every
#' step after that is arithmetic on those tables. So the vector handed to the
#' gini calculation in steps 3 and 4 holds **one value per cluster**, not one
#' per cell — its length is the number of clusters being compared. Two
#' consequences worth knowing before setting `min_expression_gini` or
#' `min_detection_gini`:
#'
#' * a gini coefficient over `G` values cannot exceed `(G - 1) / G`, so the
#'   scale depends on how many clusters are being compared: the ceiling is
#'   0.95 for 20 clusters, 0.80 for 5, and 0.50 for 2. A threshold chosen for
#'   a many-cluster run is therefore not transferable to a smaller one, where
#'   it may sit above the ceiling and reject every feature;
#' * [findGiniMarkers_one_vs_all()] compares each cluster against all others
#'   pooled into a single group, so its coefficients are always taken over
#'   exactly two values and are capped at 0.50 — well below what the same
#'   feature would score in a multi-cluster run here.
#'
#' `min_length` exists to remove that dependence; see the *Comparing runs*
#' section below.
#'
#' As a result "top gini" feats are feats that are very selectively expressed
#' in a specific cluster,
#' however not always expressed in all cells of that cluster. In other words
#' highly specific, but
#' not necessarily sensitive at the single-cell level.
#'
#' @section Filtering:
#' Two kinds of gate are available, and they do different jobs.
#'
#' `min_expression` and `min_detection` are an abundance floor. The gini
#' coefficient is scale-free: a feature averaging 0.001 in one cluster against
#' 0.0001 in the rest scores exactly as specific as one averaging 100 against
#' 10. It therefore carries no magnitude term, and will rank near-noise
#' features as perfectly selective. These two gates supply that term. The
#' equivalent for the other methods is `logFC`, which
#' [findScranMarkers_one_vs_all()] and [findMastMarkers_one_vs_all()] filter
#' on directly.
#'
#' `min_expression_gini` and `min_detection_gini` gate the gini coefficients
#' themselves, and are off by default. They are best used on a second pass:
#' the returned table always carries the `expression_gini` and `detection_gini`
#' columns, so run once, look at their distribution, then set a floor. What
#' they add over `min_feats` is an *absolute* cutoff, where `min_feats` is a
#' *relative* one -- it takes each cluster's top few by `comb_score` however
#' weak those are, while a gini floor is a statement about the coefficient
#' itself.
#'
#' `rank_score` is a third kind: a position rather than a value. It keeps a
#' feature when the cluster in question is within the given rank for both
#' `expression` and `detection`, rank 1 being the cluster where the feature is
#' highest. So `rank_score = 1` means "this cluster tops the feature on both
#' measures", and `2` means "top two". It is off by default (`Inf`). Clusters
#' tied at the top all hold the same rank, so a tie does not exclude them.
#'
#' All four value comparisons are strict (`>`), so a value exactly equal to
#' its threshold does not pass; `rank_score` is inclusive (`<=`). The strict
#' case is worth knowing for `min_detection`, where the values are fractions
#' over a cluster's cell count and land on round numbers often.
#'
#' Every gate is combined with `min_feats` by `or`, so a feature in the top
#' `min_feats` of its cluster by `comb_score` is kept regardless of all of
#' them. Setting the gates more strictly therefore shrinks the result towards
#' `min_feats` features per cluster, never below it.
#'
#' @section Comparing runs:
#' A *run* here means one call together with the set of clusters that call
#' compares — `G` of them. A gini coefficient over `G` values cannot exceed
#' `(G - 1) / G`, so a coefficient is a property of the feature **and** of the
#' run that produced it, not of the feature alone.
#'
#' `G` is not fixed by the dataset; it is set by the arguments. On the mini
#' visium object with seven leiden clusters:
#'
#' | call | `G` | ceiling | highest `expression_gini` |
#' | --- | --- | --- | --- |
#' | `findGiniMarkers()` | 7 | 0.857 | 0.484 |
#' | `findGiniMarkers(subset_clusters = 3 of them)` | 3 | 0.667 | 0.443 |
#' | `findGiniMarkers(group_1 =, group_2 =)` | 2 | 0.500 | 0.390 |
#' | `findGiniMarkers_one_vs_all()` | 2 | 0.500 | 0.309 |
#'
#' The feature `Mustn1` scores 0.484 in the full run and 0.238 in the
#' three-cluster subset — the same cells, the same expression, half the
#' coefficient, because `subset_clusters` changed what it was compared against.
#' That is the boundary: any point where a coefficient produced under one `G`
#' is read against one produced under another.
#'
#' `min_length` pads the per-cluster vector with copies of its own minimum
#' before the coefficient is taken, fixing the ceiling at
#' `(min_length - 1) / min_length` for any run with fewer clusters than that,
#' and so putting runs with different `G` back on one scale.
#'
#' Padding is not a neutral rescaling, and it is not strictly better than
#' leaving it off. Two things follow from the padding value being *each
#' feature's own minimum*. It reorders features within a single run, because a
#' feature with a floor of zero gains far more from padding than one with a
#' high floor — measured on mini visium, padded and unpadded
#' `expression_gini` correlate at 0.95, not 1.0. And it asserts something the
#' data does not contain: `c(10, 0)` scores 0.50 over two groups and 0.9375
#' padded to sixteen, which amounts to assuming the feature would sit at its
#' observed minimum in fourteen groups nobody measured. Being top of two
#' clusters really is weaker evidence of specificity than being top of twenty;
#' the unpadded ceiling is that fact, not an artefact.
#'
#' So this is a deliberate trade — a comparable number in exchange for an
#' assumption — and it cannot be chosen automatically, because the right value
#' depends on which *other* runs the score has to line up with, which a single
#' call cannot see.
#'
#' **Leave it at `0` when** you are reading one run's results on their own
#' terms — the default, and the case that needs no assumption. Nothing in the
#' returned table requires an absolute coefficient: `comb_rank`, `min_feats`
#' and `rank_score` are all relative to the run, and the gini columns are only
#' being read against each other.
#'
#' **Set it when a gini number has to mean the same thing twice.** In rough
#' order of how easily each is hit:
#'
#' * comparing a full run against one narrowed by `subset_clusters`, or against
#'   a `group_1`/`group_2` pairwise call — the same object and the same
#'   clustering, but a different `G`, so the coefficients are not on one scale;
#' * comparing [findGiniMarkers()] against [findGiniMarkers_one_vs_all()],
#'   whose two-group comparisons are capped at 0.50 while a 20-cluster run
#'   reaches 0.95;
#' * reusing one `min_expression_gini` or `min_detection_gini` threshold across
#'   clustering resolutions, or across datasets, where the cluster count
#'   differs;
#' * reporting a coefficient as a property of a feature rather than of one
#'   analysis.
#'
#' All four are the same situation: a coefficient produced under one `G` being
#' read against one produced under another.
#'
#' Pick a value at least as large as the biggest cluster count you want to
#' compare across. Below that, a run with more clusters than `min_length` is
#' left unpadded and the scales still differ. Above it there is no real cost:
#' padding harder widens the range rather than compressing it, and barely
#' touches the ordering (measured on the mini visium dataset, `min_length` of
#' 16 versus 100 gives a Spearman correlation of 0.998 between the resulting
#' coefficients). `16` is the value this padding shipped with historically and
#' fixes the ceiling at 0.9375.
#'
#' Turning it on is not free. Every gini score changes, so results stop being
#' comparable with unpadded runs — pad both sides or neither. And since
#' `comb_score` is built from the coefficients, the reordering described above
#' propagates: on mini visium `min_length = 16` moves `comb_rank` for 98% of
#' rows, which changes which features `min_feats` rescues and therefore which
#' features are returned at all.
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
        min_expression = 0.2,
        min_detection = 0.2,
        min_expression_gini = -Inf,
        min_detection_gini = -Inf,
        detection_threshold = 0,
        min_length = 0,
        rank_score = Inf,
        min_feats = 5,
        min_genes = NULL,
        min_expr_gini_score = deprecated(),
        min_det_gini_score = deprecated()) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
    }
    # `min_expr_gini_score` / `min_det_gini_score` never gated the gini
    # coefficients despite their names -- they gated `expression` and
    # `detection`, the per-cluster mean and detection fraction. The new names
    # say that, and the gini coefficients get gates of their own above.
    .dep <- function(...) {
        deprecate_param(..., fun = "findGiniMarkers", when = "4.2.4")
    }
    min_expression <- .dep(min_expr_gini_score, min_expression)
    min_detection <- .dep(min_det_gini_score, min_detection)

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


    # Narrowing is expressed as `NA` in the grouping vector rather than by
    # subsetting the object. The aggregate excludes NA-group cells on every
    # backend, so this is equivalent, avoids copying the object, and keeps a
    # disk-backed store from being materialized just to drop columns.
    grp <- as.character(cell_metadata[][[cluster_column]])

    if (!is.null(subset_clusters)) {
        grp[!grp %in% as.character(subset_clusters)] <- NA
    } else if (!is.null(group_1) & !is.null(group_2)) {
        if (!is.null(group_1_name)) {
            if (!is.character(group_1_name)) {
                stop("group_1_name needs to be a character")
            }
        } else {
            group_1_name <- paste0(group_1, collapse = "_")
        }

        if (!is.null(group_2_name)) {
            if (!is.character(group_2_name)) {
                stop("group_2_name needs to be a character")
            }
        } else {
            group_2_name <- paste0(group_2, collapse = "_")
        }

        grp <- ifelse(grp %in% as.character(group_1), group_1_name,
            ifelse(grp %in% as.character(group_2), group_2_name, NA_character_)
        )
    }

    if (all(is.na(grp))) {
        stop("no cells remain after the requested cluster selection")
    }

    # Dispatched on the expression object, not computed here: `getExpression`
    # returns whatever the slot holds, so a disk-backed store arrives intact
    # and its own featStats method supplies the pass. Everything from the
    # statistic onward lives in analyzeData(x, giniMarkersParam).
    #
    # Cluster labels come straight from the grouping vector, so the 'cluster_'
    # prefix that create_average_DT() added -- and the conditional strip that
    # undid it -- are both gone. That strip only fired when no real cluster name
    # contained "cluster_", so a cluster genuinely named 'cluster_3' used to
    # come back mangled; labels are now verbatim.
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )
    analyzeData(
        x = expr_data,
        param = markersParam(
            method = "gini",
            comparison = "pairwise",
            min_expression = min_expression,
            min_detection = min_detection,
            min_expression_gini = min_expression_gini,
            min_detection_gini = min_detection_gini,
            detection_threshold = detection_threshold,
            min_length = min_length,
            rank_score = rank_score,
            min_feats = min_feats
        ),
        # named so the verb matches on cell ID rather than position
        groups = stats::setNames(grp, cell_metadata[][["cell_ID"]])
    )
}


# analyzeData(<any>, giniMarkersParam) ####

# Gini marker detection consumes `analyzeData(x, featStatsParam)` and nothing
# else, so it is substrate-agnostic by construction -- hence the `ANY`
# signature. There is deliberately no per-backend method: a store that
# implements the grouped statistic supplies gini markers already, which is the
# whole reason the statistic was pushed into that verb.
#
# `ANY` is the primary implementation here, not a warning fallback as it is for
# `autoPcaParam`. What `x` must satisfy is checked below rather than encoded in
# the signature, because the requirement is "featStats dispatches on it", which
# a class union cannot express across package boundaries.

#' @rdname markers_gini
#' @param x expression values — anything `analyzeData(x, featStatsParam)`
#'   accepts, including a `matrix`, a `Matrix`, a `DelayedMatrix` (what
#'   `expression_values = "scaled"` holds) or a disk-backed store.
#' @param param a [giniMarkersParam-class].
#' @param groups vector of group assignments, one per cell. Taken in column
#'   order of `x`, unless the vector is named by cell ID, in which case it is
#'   matched on identity. Naming is the safer form: the expression matrix and
#'   the cell metadata a caller builds `groups` from are fetched independently
#'   and need not share a cell order. `NA` excludes a cell from every group,
#'   which is how a caller narrows to a subset without copying the object.
#' @param verbose report progress per group. `"one_vs_rest"` only.
#' @export
setMethod("analyzeData",
    signature(x = "ANY", param = "giniMarkersParam"),
    function(x, param, ..., groups = NULL, verbose = TRUE) {
        .markers_gini(x, param, groups = groups, verbose = verbose)
    }
)


#' @keywords internal
#' @noRd
.markers_gini <- function(x, param, groups, verbose = TRUE) {
    if (is.null(groups)) {
        stop("[analyzeData(giniMarkersParam)] `groups` is required: one group ",
            "assignment per cell.", call. = FALSE)
    }
    if (!hasMethod("analyzeData", signature(class(x)[1L], "featStatsParam"))) {
        stop("[analyzeData(giniMarkersParam)] no ",
            "analyzeData(<", class(x)[1L], ">, featStatsParam) method, so the ",
            "per-group statistics gini needs cannot be computed. Gini markers ",
            "are derived entirely from that verb.",
            call. = FALSE
        )
    }

    # The one pass over the values. `mean_expr` is the per-group mean and
    # `perc_cells / 100` the detection fraction; `sum` and `nnz` are the only
    # accumulators either needs, so the other two are never computed.
    st <- analyzeData(
        x,
        analyzeParam("feat_stats",
            detection_threshold = param$detection_threshold),
        groups = groups,
        stats = c("sum", "nnz")
    )

    if (identical(param$comparison %null% "pairwise", "one_vs_rest")) {
        return(.markers_one_vs_rest_gini(st, param, verbose = verbose))
    }
    .gini_score_dt(
        data.table::data.table(
            feats = st$feats,
            cluster = st$group,
            expression = st$mean_expr,
            detection = st$perc_cells / 100
        ),
        min_length = param$min_length,
        min_expression = param$min_expression,
        min_detection = param$min_detection,
        min_expression_gini = param$min_expression_gini,
        min_detection_gini = param$min_detection_gini,
        rank_score = param$rank_score,
        min_feats = param$min_feats
    )
}


# One table per group, each scoring that group against the pooled remainder.
#
# Takes the grouped statistics rather than the store: the accumulators behind
# `total_expr` and `nr_cells` are additive, so each group-vs-rest pair is a
# row-sum over the other columns. What used to be one full scan per group is one
# scan total, and the pooling is exact rather than an approximation.
#' @keywords internal
#' @noRd
.markers_one_vs_rest_gini <- function(st, param, verbose = TRUE) {
    cluster <- NULL   # NSE binding

    lvls <- unique(st$group)
    n_feats <- length(unique(st$feats))
    feat_ids <- st$feats[seq_len(n_feats)]

    # groups vary slowest with feats cycling within, so these unroll directly
    sums <- matrix(st$total_expr, nrow = n_feats, dimnames = list(NULL, lvls))
    nnz <- matrix(
        as.numeric(st$nr_cells), nrow = n_feats, dimnames = list(NULL, lvls)
    )
    n_k <- st$n_cells[seq(1L, by = n_feats, length.out = length(lvls))]
    names(n_k) <- lvls

    uniq_clusters <- mixedsort(lvls)

    with_pbar({
        pb <- pbar(along = uniq_clusters)
        result_list <- lapply(
            seq_along(uniq_clusters),
            function(clus_i) {
                selected_clus <- as.character(uniq_clusters[clus_i])
                other_clus <- setdiff(lvls, selected_clus)

                if (isTRUE(verbose)) {
                    cat("start with cluster ", selected_clus)
                }

                # default group naming matches what findGiniMarkers() would
                # have produced for group_1 = selected, group_2 = rest
                rest_name <- paste0(
                    mixedsort(uniq_clusters[
                        as.character(uniq_clusters) != selected_clus
                    ]),
                    collapse = "_"
                )
                n_rest <- sum(n_k[other_clus])

                markers <- .gini_score_dt(
                    data.table::data.table(
                        feats = rep(feat_ids, 2L),
                        cluster = rep(
                            c(selected_clus, rest_name),
                            each = n_feats
                        ),
                        expression = c(
                            sums[, selected_clus] / n_k[[selected_clus]],
                            rowSums(sums[, other_clus, drop = FALSE]) / n_rest
                        ),
                        detection = c(
                            nnz[, selected_clus] / n_k[[selected_clus]],
                            rowSums(nnz[, other_clus, drop = FALSE]) / n_rest
                        )
                    ),
                    min_length = param$min_length,
                    min_expression = param$min_expression,
                    min_detection = param$min_detection,
                    min_expression_gini = param$min_expression_gini,
                    min_detection_gini = param$min_detection_gini,
                    rank_score = param$rank_score,
                    min_feats = param$min_feats
                )

                filtered_table <- markers[cluster == selected_clus]

                pb(message = c("cluster ", clus_i, "/", length(uniq_clusters)))
                return(filtered_table)
            }
        )
    })

    do.call("rbind", result_list)
}


# Score a features x clusters table into the gini marker result.
#
# Takes `feats`, `cluster`, `expression` and `detection` and returns the
# filtered, ordered result. Shared by the pairwise and pooled one-vs-rest paths
# so the statistic has one implementation.
.gini_score_dt <- function(aggr_sc,
    min_length, min_expression, min_detection,
    min_expression_gini, min_detection_gini, rank_score, min_feats) {
    # data.table variables
    feats <- cluster <- expression <- detection <- NULL
    expression_gini <- detection_gini <- NULL
    expression_rank <- detection_rank <- expression_wt <- detection_wt <- NULL
    comb_score <- comb_rank <- NULL

    # `min_length` pads the per-cluster vector so the coefficient stops
    # depending on how many clusters were compared -- see mygini_fun(). 0, the
    # default, never pads.
    aggr_sc[, expression_gini := mygini_fun(
        expression,
        min_length = min_length
    ), by = feats]
    aggr_sc[, detection_gini := mygini_fun(
        detection,
        min_length = min_length
    ), by = feats]

    # Two distinct things are wanted from the same ordering, so they get two
    # columns. The rank is what `rank_score` filters on -- position 1 is the
    # cluster where this feat is highest. The weight is that rank rescaled
    # within the cluster to [1, 0.1], and only feeds `comb_score`. Collapsing
    # them, as this did previously, left `rank_score` comparing against a
    # value that can never exceed 1, so it never had any effect at its own
    # default.
    #
    # `ties.method = "min"` on the rank so that clusters tied at the top all
    # hold position 1. The default "average" gives every one of them 1.5 and
    # `<= 1` then rejects the feat outright rather than crediting each winner.
    # Detection is a fraction over a cluster's cell count, so ties are common.
    # The weight keeps the default tie handling, which is what `comb_score`
    # has always been built from.
    aggr_sc[, expression_rank := rank(-expression, ties.method = "min"),
        by = feats]
    aggr_sc[, expression_wt := rank(-expression), by = feats]
    aggr_sc[, expression_wt := scales::rescale(
        expression_wt,
        to = c(1, 0.1)
    ), by = cluster]

    aggr_sc[, detection_rank := rank(-detection, ties.method = "min"),
        by = feats]
    aggr_sc[, detection_wt := rank(-detection), by = feats]
    aggr_sc[, detection_wt := scales::rescale(
        detection_wt,
        to = c(1, 0.1)
    ), by = cluster]

    aggr_sc[, comb_score := (expression_gini * expression_wt) * (
        detection_gini * detection_wt)]
    setorder(aggr_sc, cluster, -comb_score)
    aggr_sc[, comb_rank := seq_len(.N), by = cluster]

    top_feats_scores <- aggr_sc[comb_rank <= min_feats | (
        expression_rank <= rank_score & detection_rank <= rank_score)]
    # Gini is scale-free -- a 0.001 vs 0.0001 difference between clusters
    # scores identically to 100 vs 10 -- so it carries no magnitude term of its
    # own and will rank near-noise features as perfectly specific. The
    # expression and detection floors are that magnitude term; scran and MAST
    # get theirs for free from logFC. The gini floors are the absolute cutoff
    # `comb_rank` cannot express, since it only ever ranks within a cluster.
    top_feats_scores_filtered <- top_feats_scores[comb_rank <= min_feats | (
        expression > min_expression &
            detection > min_detection &
            expression_gini > min_expression_gini &
            detection_gini > min_detection_gini)]
    setorder(top_feats_scores_filtered, cluster, comb_rank)

    # the rescaled weights exist only to build `comb_score`, which is returned;
    # drop them so the column contract is unchanged
    top_feats_scores_filtered[, c("expression_wt", "detection_wt") := NULL]

    # stated rather than left to fall out of construction order
    data.table::setcolorder(top_feats_scores_filtered, c(
        "feats", "cluster", "expression", "expression_gini",
        "detection", "detection_gini", "expression_rank",
        "detection_rank", "comb_score", "comb_rank"
    ))

    top_feats_scores_filtered[]
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
#' @param min_expression minimum per-cluster mean expression, gating the
#' `expression` column of the result
#' @param min_detection minimum fraction of a cluster's cells with expression
#' above `detection_threshold`, gating the `detection` column of the result
#' @param min_expression_gini minimum gini coefficient of expression, gating
#' the `expression_gini` column of the result. `-Inf` (default) disables it.
#' @param min_detection_gini minimum gini coefficient of detection, gating the
#' `detection_gini` column of the result. `-Inf` (default) disables it.
#' @param detection_threshold expression value above which a cell counts as
#' expressing a feature, used when computing `detection`. Not a filter on the
#' returned rows -- see `min_detection` for that.
#' @param min_length pad the per-cluster vector to this length before taking
#' the gini coefficient, using copies of its minimum. Removes the dependence
#' of the coefficient on how many clusters were compared, so gini scores and
#' the `min_expression_gini` / `min_detection_gini` thresholds become
#' comparable across runs. `0` (the default) never pads.
#' @param rank_score keep a feature when its cluster is within this rank for
#' both `expression` and `detection`, where rank 1 is the cluster in which the
#' feature is highest. `Inf` (default) disables it. Combined with `min_feats`
#' by `or`, like the other gates.
#' @param min_feats minimum number of top feats to return
#' @param min_genes deprecated, use min_feats
#' @param verbose be verbose
#' @param min_expr_gini_score `r lifecycle::badge("deprecated")` use
#' `min_expression`. Despite its name it never gated a gini coefficient.
#' @param min_det_gini_score `r lifecycle::badge("deprecated")` use
#' `min_detection`. Despite its name it never gated a gini coefficient.
#' @returns data.table with marker feats
#' @details
#' Each cluster is compared against every other cluster pooled into a single
#' group, by calling [findGiniMarkers()] once per cluster and keeping the rows
#' belonging to the cluster under test. See there for how the scores are built.
#'
#' Because each of those calls compares exactly two groups, the gini
#' coefficients are taken over two values and so cannot exceed 0.50 — the
#' ceiling for `G` groups is `(G - 1) / G`. Thresholds passed to
#' `min_expression_gini` or `min_detection_gini` have to sit below that or they
#' reject every feature, leaving only the `min_feats` per cluster that the
#' filter always keeps.
#'
#' Two groups also gives `rank_score` a direct reading here: with only the
#' cluster and the pooled rest to rank, `rank_score = 1` keeps a feature only
#' where the cluster under test is not beaten by the rest on either mean
#' expression or detection fraction — a tie counts, since tied groups share
#' rank 1. Any value above 1 admits every feature, there being no third rank
#' to exclude.
#' @md
#' @inheritSection findGiniMarkers Filtering
#' @inheritSection findGiniMarkers Comparing runs
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
        min_expression = 0.5,
        min_detection = 0.5,
        min_expression_gini = -Inf,
        min_detection_gini = -Inf,
        detection_threshold = 0,
        min_length = 0,
        rank_score = Inf,
        min_feats = 4,
        min_genes = NULL,
        verbose = TRUE,
        min_expr_gini_score = deprecated(),
        min_det_gini_score = deprecated()) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
    }
    .dep <- function(...) {
        deprecate_param(
            ..., fun = "findGiniMarkers_one_vs_all", when = "4.2.4"
        )
    }
    min_expression <- .dep(min_expr_gini_score, min_expression)
    min_detection <- .dep(min_det_gini_score, min_detection)

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


    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )
    analyzeData(
        x = expr_data,
        param = markersParam(
            method = "gini",
            comparison = "one_vs_rest",
            min_expression = min_expression,
            min_detection = min_detection,
            min_expression_gini = min_expression_gini,
            min_detection_gini = min_detection_gini,
            detection_threshold = detection_threshold,
            min_length = min_length,
            rank_score = rank_score,
            min_feats = min_feats
        ),
        # named so the verb matches on cell ID rather than position
        groups = stats::setNames(
            as.character(cell_metadata[[cluster_column]]),
            cell_metadata[["cell_ID"]]
        ),
        verbose = verbose
    )
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
#' @param min_expression gini: minimum per-cluster mean expression
#' @param min_detection gini: minimum fraction of a cluster's cells with
#' expression above `detection_threshold`
#' @param min_expression_gini gini: minimum gini coefficient of expression.
#' `-Inf` (default) disables it.
#' @param min_detection_gini gini: minimum gini coefficient of detection.
#' `-Inf` (default) disables it.
#' @param detection_threshold gini: expression value above which a cell counts
#' as expressing a feature
#' @param min_length gini: pad the per-cluster vector to this length before
#' taking the gini coefficient, making scores comparable across runs with
#' different cluster counts. `0` (the default) never pads.
#' @param rank_score gini: keep a feature when its cluster is within this
#' rank for both `expression` and `detection` (rank 1 = highest). `Inf`
#' (default) disables it.
#' @param min_feats minimum number of top feats to return (for gini)
#' @param min_genes deprecated, use min_feats
#' @param group_1_name mast: custom name for group_1 clusters
#' @param group_2_name mast: custom name for group_2 clusters
#' @param adjust_columns mast: column in pDataDT to adjust for
#' (e.g. detection rate)
#' @param min_expr_gini_score `r lifecycle::badge("deprecated")` use
#' `min_expression`. Despite its name it never gated a gini coefficient.
#' @param min_det_gini_score `r lifecycle::badge("deprecated")` use
#' `min_detection`. Despite its name it never gated a gini coefficient.
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

# TODO: findMarkers() and analyzeData(x, markersParam) occupy the same mental
# space, and keeping both means every marker method has two entry points with
# different argument vocabularies and different defaults (see the 0.2/0.2 vs
# 0.5/0.5 split between findGiniMarkers() and this function). Pick one:
# either deprecate findMarkers() in favour of param dispatch, or promote it to
# a generic of its own rather than layering it over analyzeData(). Doing
# neither is what lets the two drift.
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
        min_expression = 0.5,
        min_detection = 0.5,
        min_expression_gini = -Inf,
        min_detection_gini = -Inf,
        detection_threshold = 0,
        min_length = 0,
        rank_score = Inf,
        min_feats = 4,
        min_genes = NULL,
        group_1_name = NULL,
        group_2_name = NULL,
        adjust_columns = NULL,
        min_expr_gini_score = deprecated(),
        min_det_gini_score = deprecated(),
        ...) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
    }
    .dep <- function(...) {
        deprecate_param(..., fun = "findMarkers", when = "4.2.4")
    }
    min_expression <- .dep(min_expr_gini_score, min_expression)
    min_detection <- .dep(min_det_gini_score, min_detection)

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
            min_expression = min_expression,
            min_detection = min_detection,
            min_expression_gini = min_expression_gini,
            min_detection_gini = min_detection_gini,
            detection_threshold = detection_threshold,
            min_length = min_length,
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
#' @param min_expression gini: minimum per-cluster mean expression
#' @param min_detection gini: minimum fraction of a cluster's cells with
#' expression above `detection_threshold`
#' @param min_expression_gini gini: minimum gini coefficient of expression.
#' `-Inf` (default) disables it.
#' @param min_detection_gini gini: minimum gini coefficient of detection.
#' `-Inf` (default) disables it.
#' @param detection_threshold gini: expression value above which a cell counts
#' as expressing a feature
#' @param min_length gini: pad the per-cluster vector to this length before
#' taking the gini coefficient, making scores comparable across runs with
#' different cluster counts. `0` (the default) never pads.
#' @param rank_score gini: keep a feature when its cluster is within this
#' rank for both `expression` and `detection` (rank 1 = highest). `Inf`
#' (default) disables it.
#' @param adjust_columns mast: column in pDataDT to adjust for
#' (e.g. detection rate)
#' @param verbose be verbose
#' @param min_expr_gini_score `r lifecycle::badge("deprecated")` use
#' `min_expression`. Despite its name it never gated a gini coefficient.
#' @param min_det_gini_score `r lifecycle::badge("deprecated")` use
#' `min_detection`. Despite its name it never gated a gini coefficient.
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
        min_expression = 0.5,
        min_detection = 0.5,
        min_expression_gini = -Inf,
        min_detection_gini = -Inf,
        detection_threshold = 0,
        min_length = 0,
        rank_score = Inf,
        # mast specific
        adjust_columns = NULL,
        verbose = TRUE,
        min_expr_gini_score = deprecated(),
        min_det_gini_score = deprecated(),
        ...) {
    ## deprecated arguments
    if (!is.null(min_genes)) {
        min_feats <- min_genes
        warning("min_genes argument is deprecated, use min_feats argument in
                the future")
    }
    .dep <- function(...) {
        deprecate_param(..., fun = "findMarkers_one_vs_all", when = "4.2.4")
    }
    min_expression <- .dep(min_expr_gini_score, min_expression)
    min_detection <- .dep(min_det_gini_score, min_detection)

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
            min_expression = min_expression,
            min_detection = min_detection,
            min_expression_gini = min_expression_gini,
            min_detection_gini = min_detection_gini,
            detection_threshold = detection_threshold,
            min_length = min_length,
            rank_score = rank_score,
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
