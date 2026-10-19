# ============================================================================
# Spatial co-expression modules.
#
# Genes whose spatially smoothed expression is correlated, grouped into
# modules. spatialCoexpression() composes generic building blocks, each of
# which dispatches on the expression object, so a disk-backed store works
# wherever its backend provides them:
#
#   processData(x, spatialSmoothParam, network =)  # GiottoClass
#   analyzeData(x, featCorParam)                   # GiottoClass
#   hclust on 1 - correlation                      # modules
#   createMetafeats()                              # one score per module
# ============================================================================


#' @title Spatial co-expression modules
#' @name spatialCoexpression
#' @description
#' Groups features into modules whose expression varies together across the
#' tissue. Each cell's values are first averaged with its neighbours' in a
#' spatial network ([GiottoClass::spatialSmoothParam()]), the features are
#' correlated on the smoothed values, and the correlation is clustered
#' hierarchically into `n_modules` modules. Each module then gets one score
#' per cell, the mean of its features ([GiottoClass::createMetafeats()]).
#' @param gobject giotto object
#' @param spat_unit,feat_type spatial unit and feature type
#' @param expression_values expression values to smooth and score
#' @param spatial_network_name spatial network to smooth over
#' @param feats features to group, typically spatially variable features from
#' [spatialAutoCor()]. Default all features; the correlation matrix is
#' features x features, so keep this to a few thousand at most.
#' @param smoothing weight of the neighbours in `[0, 1]`; `NULL` uses
#' \eqn{1 - 1/k} with \eqn{k} the median number of neighbours, and `0` uses
#' the unsmoothed values (plain co-expression)
#' @param weighted weight neighbours by the network's `weight` column
#' @param hops average over neighbours up to this many edges away
#' @param cor_method `"pearson"` or `"spearman"`
#' @param n_modules number of modules
#' @param cluster_method agglomeration method passed to [stats::hclust()]
#' @param score_stat statistic for the module scores, see
#' [GiottoClass::createMetafeats()]; `"zmean"` keeps highly expressed
#' features from dominating
#' @param name name of the feature metadata column with the module of each
#' feature, and of the module scores
#' @param return_gobject return the giotto object with the results added
#' (default), or the results as a list
#' @param verbose be verbose
#' @returns With `return_gobject = TRUE`, the giotto object with, in the
#' feature metadata, `<name>` (module of each feature, `NA` for features not
#' grouped) and `<name>_cor` (the feature's mean correlation with the rest of
#' its module), and module scores per cell as the spatial enrichment
#' `<name>`. Otherwise a list:
#' * `cor`: the feature x feature correlation of the smoothed values
#' * `modules`: module of each feature (named integer vector)
#' * `feats`: data.table of `feat_ID`, `module` and `module_cor`
#' * `summary`: one row per module with its size, `within_cor` (mean
#'   correlation between its features), `between_cor` (mean correlation with
#'   features of other modules) and `score`, which ranks modules that are
#'   coherent and distinct from the others first
#' * `hclust`: the clustering
#'
#' Features with the same value in every cell have no correlation and are
#' left out of the modules.
#' @seealso [heatmapCoexpression()]
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' sv <- spatialAutoCor(g, method = "binspect",
#'     spatial_network_name = "spatial_network", return_gobject = FALSE)
#' g <- spatialCoexpression(g, feats = sv$feats[1:100],
#'     spatial_network_name = "spatial_network", n_modules = 5)
#' spatCellPlot(g, spat_enr_names = "spatial_module",
#'     cell_annotation_values = as.character(1:5))
#' @export
spatialCoexpression <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = "normalized",
    spatial_network_name = "Delaunay_network",
    feats = NULL,
    smoothing = NULL,
    weighted = FALSE,
    hops = 1L,
    cor_method = c("pearson", "spearman"),
    n_modules = 10L,
    cluster_method = "ward.D",
    score_stat = c("mean", "zmean"),
    name = "spatial_module",
    return_gobject = TRUE,
    verbose = NULL) {
    i.module <- i.module_cor <- NULL  # NSE
    cor_method <- match.arg(cor_method)
    score_stat <- match.arg(score_stat)
    checkmate::assert_count(n_modules, positive = TRUE)
    spat_unit <- set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject = gobject,
        spat_unit = spat_unit, feat_type = feat_type)

    x <- getExpression(gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = expression_values, output = "exprObj")[]
    if (!is.null(feats)) {
        miss <- setdiff(feats, rownames(x))
        if (length(miss)) {
            stop("[spatialCoexpression] features not found: ",
                paste(head(miss, 5), collapse = ", "), call. = FALSE)
        }
        x <- x[feats, , drop = FALSE]
    }
    net <- getSpatialNetwork(gobject, spat_unit = spat_unit,
        name = spatial_network_name, output = "spatialNetworkObj",
        verbose = FALSE)[]

    vmsg(.v = verbose, "smoothing", nrow(x), "features over",
        spatial_network_name)
    sm <- processData(x, spatialSmoothParam(smoothing = smoothing,
        weighted = weighted, hops = hops), network = net)
    cm <- analyzeData(sm, featCorParam(method = cor_method))
    res <- .coexpr_modules(cm, n_modules, cluster_method)
    if (!isTRUE(return_gobject)) return(res)

    fm <- fDataDT(gobject, spat_unit = spat_unit,
        feat_type = feat_type)[, "feat_ID"]
    fm[res$feats, on = "feat_ID", (name) := i.module]
    fm[res$feats, on = "feat_ID", (paste0(name, "_cor")) := i.module_cor]
    gobject <- addFeatMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, new_metadata = fm, by_column = TRUE,
        column_feat_ID = "feat_ID")
    createMetafeats(gobject, spat_unit = spat_unit, feat_type = feat_type,
        expression_values = expression_values,
        feat_clusters = res$modules, stat = score_stat, name = name,
        verbose = FALSE)
}


#' @title Heatmap of spatial co-expression modules
#' @name heatmapCoexpression
#' @description Correlation heatmap of the features grouped by
#' [spatialCoexpression()], ordered by its clustering and annotated by
#' module.
#' @param gobject giotto object, for the plot instructions
#' @param coexpression the list returned by
#' `spatialCoexpression(return_gobject = FALSE)`
#' @param show_row_names,show_column_names show feature names
#' @param show_row_dend show the row dendrogram
#' @inheritParams GiottoVisuals::plot_output_params
#' @param ... additional arguments to [ComplexHeatmap::Heatmap()]
#' @returns a ComplexHeatmap heatmap, depending on the plot settings
#' @export
heatmapCoexpression <- function(gobject,
    coexpression,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = TRUE,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "heatmapCoexpression",
    ...) {
    package_check(pkg_name = "ComplexHeatmap", repository = "CRAN")
    if (!all(c("cor", "modules", "hclust") %in% names(coexpression))) {
        stop("[heatmapCoexpression] `coexpression` must be the result of ",
            "spatialCoexpression(return_gobject = FALSE)", call. = FALSE)
    }
    keep <- names(coexpression$modules)
    m <- coexpression$cor[keep, keep, drop = FALSE]
    mods <- as.character(coexpression$modules)
    cols <- getDistinctColors(length(unique(mods)))
    names(cols) <- unique(mods)
    ha <- ComplexHeatmap::HeatmapAnnotation(module = mods,
        col = list(module = cols))
    heatm <- ComplexHeatmap::Heatmap(m,
        name = "correlation",
        cluster_rows = coexpression$hclust,
        cluster_columns = coexpression$hclust,
        show_row_dend = show_row_dend, show_column_dend = FALSE,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        top_annotation = ha, ...)
    GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = heatm,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    )
}


# internals ####

# Modules from a feature x feature correlation: hclust on 1 - r, cut into
# `k`. Features without a defined correlation are left out.
.coexpr_modules <- function(cm, k, method) {
    score <- within_cor <- between_cor <- NULL  # NSE
    ok <- !is.na(diag(cm))
    r <- cm[ok, ok, drop = FALSE]
    r[is.na(r)] <- 0
    k <- min(k, nrow(r))
    hc <- stats::hclust(stats::as.dist(1 - r), method = method)
    modules <- stats::cutree(hc, k = k)

    # each feature's mean correlation with the rest of its module
    module_cor <- vapply(seq_along(modules), function(i) {
        same <- setdiff(which(modules == modules[i]), i)
        if (length(same)) mean(r[i, same]) else NA_real_
    }, numeric(1))
    feats <- data.table::data.table(feat_ID = names(modules),
        module = unname(modules), module_cor = module_cor)

    summary <- data.table::rbindlist(lapply(sort(unique(modules)), function(m) {
        i <- modules == m
        w <- r[i, i, drop = FALSE]
        data.table::data.table(module = m, size = sum(i),
            within_cor = if (sum(i) > 1) mean(w[upper.tri(w)]) else NA_real_,
            between_cor = if (any(!i)) mean(r[i, !i]) else NA_real_)
    }))
    # coherent modules that anticorrelate with the rest rank first
    summary[, score := within_cor *
        (1 - (between_cor - min(between_cor, na.rm = TRUE)))]
    data.table::setorder(summary, -score, na.last = TRUE)

    # `[]` so the first print after `:=` is not swallowed
    list(cor = cm, modules = modules, feats = feats[], summary = summary[],
        hclust = hc)
}
