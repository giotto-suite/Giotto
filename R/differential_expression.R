# TODO: remove both helpers below once analyzeData(x, featStatsParam) reaches
# this branch from gsource. They are a temporary local copy: the grouped verb
# computes both statistics in one pass -- `mean_expr` is the average and
# `perc_cells / 100` the detection fraction -- and works on disk-backed
# expression, which these do not.
#
# GiottoClass exported equivalents of these. `create_average_detection_DT()`
# has since been removed there, and `create_average_DT()` is retained only for
# an internal GiottoClass caller, so Giotto now carries its own.
#
# Group membership is keyed on `cell_ID`: `getExpression()` and
# `getCellMetadata()` are fetched independently and the suite makes no
# guarantee that they share a cell order.
.average_by_group <- function(gobject, spat_unit, feat_type, meta_data_name,
    expression_values, detection_threshold = NULL) {
    expr_data <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = expression_values,
        output = "matrix"
    )
    cell_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )
    if (!meta_data_name %in% colnames(cell_metadata)) {
        stop("metadata column '", meta_data_name, "' not found",
            call. = FALSE
        )
    }

    ord <- match(colnames(expr_data), cell_metadata[["cell_ID"]])
    if (anyNA(ord)) {
        stop(
            "expression columns and cell metadata do not describe the same ",
            "cells; cannot align them.",
            call. = FALSE
        )
    }
    grouping <- cell_metadata[[meta_data_name]][ord]

    detection <- !is.null(detection_threshold)
    savelist <- lapply(unique(cell_metadata[[meta_data_name]]), function(group) {
        temp <- expr_data[, grouping == group, drop = FALSE]
        if (detection) {
            rowSums_flex(as.matrix(temp) > detection_threshold) / ncol(temp)
        } else {
            rowMeans_flex(temp)
        }
    })
    names(savelist) <- paste0(
        "cluster_", unique(cell_metadata[[meta_data_name]])
    )

    out <- do.call("cbind", savelist)
    rownames(out) <- rownames(expr_data)
    as.data.frame(out)
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
    aggr_sc_clusters <- .average_by_group(
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
    aggr_detection_sc_clusters <- .average_by_group(
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

    # `min_length` pads the per-cluster vector so the coefficient stops
    # depending on how many clusters were compared -- see mygini_fun(). 0, the
    # default, never pads.
    aggr_sc_clusters_DT_melt[, expression_gini := mygini_fun(
        expression,
        min_length = min_length
    ), by = feats]
    aggr_detection_sc_clusters_DT_melt[, detection_gini := mygini_fun(
        detection,
        min_length = min_length
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
    expression_wt <- detection_wt <- NULL

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

    # detection rank for each feat in all samples
    # rescale detection rank range between 1 and 0.1
    aggr_sc[, detection_rank := rank(-detection, ties.method = "min"),
        by = feats]
    aggr_sc[, detection_wt := rank(-detection), by = feats]
    aggr_sc[, detection_wt := scales::rescale(
        detection_wt,
        to = c(1, 0.1)
    ), by = cluster]

    # create combine score based on rescaled ranks and gini scores

    # data.table variables
    comb_score <- comb_rank <- NULL

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
                    min_expression = min_expression,
                    min_detection = min_detection,
                    min_expression_gini = min_expression_gini,
                    min_detection_gini = min_detection_gini,
                    detection_threshold = detection_threshold,
                    min_length = min_length,
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
