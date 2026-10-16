#' @include spatautocor-param.R
NULL

# ============================================================================
# Spatial autocorrelation as compositions of building blocks:
#
#   [binSpect] processData(x, thresholdParam("binarize", <rule>))
#              analyzeData(x, edgeStatsParam(), network =)   (GiottoClass)
#              .moran_test / .geary_test / .joincount_test   (spatautocor-stats.R)
#
# Dispatching on ANY means whatever backend implements the two blocks -- a
# matrix, a GiottoDisk store -- gets every statistic, with no method of its
# own. The giotto method only resolves inputs and writes results back.
# ============================================================================

#' @rdname spatAutoCorParam
#' @param x data to test: a feature x cell matrix, an expression store, or a
#'   `giotto` object
#' @param param a [spatAutoCorParam-class]
#' @param network the spatial network matching `x`: for a matrix, a
#'   `data.frame` of `from`, `to` (cell IDs) and optionally `weight`, or an
#'   `igraph`; for a GiottoDisk store, its `parquetEdgeStore`
#' @returns `analyzeData()`: a `data.table` with one row per feature, ordered
#'   by decreasing `z`: `feats`, `statistic`, `expected`, `variance`, `z`,
#'   `p.value`, `adj.p.value`; binSpect adds `n_on` (cells called on) and
#'   `odds_ratio` (of the 2x2 table of edge types). With
#'   `return_gobject = TRUE`, the `giotto` object with `<method>.z`,
#'   `<method>.p` and `<method>.padj` added to the feature metadata, plus
#'   `binSpect.pval` (the adjusted p-value) for binSpect.
#' @export
setMethod("analyzeData",
    signature(x = "ANY", param = "spatAutoCorParam"),
    function(x, param, ..., network) {
        if (is(param, "binSpectParam")) {
            x <- processData(x, thresholdParam("binarize",
                threshold = param$bin_method,
                percentage_rank = param$percentage_rank))
        }
        st <- analyzeData(x, edgeStatsParam(weighted = isTRUE(param$weighted)),
            network = network)
        .spat_autocor_finish(param, st)
    }
)

#' @rdname spatAutoCorParam
#' @param spat_unit,feat_type spatial unit and feature type
#' @param expression_values expression values to use
#' @param spatial_network_name name of the spatial network to use
#' @param return_gobject logical. Return the `giotto` object with the results
#'   in its feature metadata instead of the result table.
#' @export
setMethod("analyzeData",
    signature(x = "giotto", param = "spatAutoCorParam"),
    function(x, param, ...,
    spat_unit = NULL, feat_type = NULL,
    expression_values = "normalized",
    spatial_network_name = "Delaunay_network",
    return_gobject = FALSE) {
        spat_unit <- set_default_spat_unit(gobject = x, spat_unit = spat_unit)
        feat_type <- set_default_feat_type(
            gobject = x, spat_unit = spat_unit, feat_type = feat_type
        )
        expr <- getExpression(x,
            values = expression_values, spat_unit = spat_unit,
            feat_type = feat_type, output = "exprObj"
        )[]
        net <- getSpatialNetwork(x,
            spat_unit = spat_unit, name = spatial_network_name,
            output = "spatialNetworkObj", verbose = FALSE
        )[]
        res <- analyzeData(expr, param, network = net)
        if (!isTRUE(return_gobject)) return(res)

        prefix <- param$method
        meta <- data.table::data.table(feats = res$feats, z = res$z,
            p = res$p.value, padj = res$adj.p.value)
        data.table::setnames(meta, c("z", "p", "padj"),
            paste0(prefix, c(".z", ".p", ".padj")))
        if (is(param, "binSpectParam")) meta[, "binSpect.pval" := res$adj.p.value]
        addFeatMetadata(x,
            spat_unit = spat_unit, feat_type = feat_type,
            new_metadata = meta, by_column = TRUE, column_feat_ID = "feats"
        )
    }
)

# statistic for the param class, adjusted p-values, ordered by evidence
.spat_autocor_finish <- function(param, st) {
    z <- p.value <- NULL # NSE
    res <- if (is(param, "moranParam")) {
        .moran_test(st)
    } else if (is(param, "gearyParam")) {
        .geary_test(st)
    } else {
        .joincount_test(st)
    }
    res[, "adj.p.value" := stats::p.adjust(p.value,
        method = param$adjust_method %null% "fdr")]
    data.table::setorder(res, -z, na.last = TRUE)
    res[]
}

# Moran's I or Geary's C with analytic significance for a cells x features
# value table and an N x N spatial weight matrix, as spatialAutoCorGlobal()
# supplies them. Both statistics and their randomisation moments are unchanged
# by replacing W with its symmetric part (W + t(W)) / 2, which is what an
# undirected weighted edge list stands for, so the weight matrix becomes one
# and goes through the same blocks as everything else.
.spat_autocor_weights <- function(use_values, feats, weight_matrix, method) {
    vals <- if (inherits(use_values, "data.table")) {
        as.matrix(use_values[, feats, with = FALSE])
    } else {
        as.matrix(use_values[, feats, drop = FALSE])
    }
    N <- nrow(vals)
    ids <- as.character(seq_len(N))
    X <- Matrix::Matrix(t(vals), sparse = TRUE)
    dimnames(X) <- list(feats, ids)

    W <- methods::as(weight_matrix, "CsparseMatrix")
    if (any(Matrix::diag(W) != 0)) {
        warning("[spatial autocorrelation] the weight matrix has a nonzero ",
            "diagonal; self-weights are ignored.", call. = FALSE)
    }
    Ws <- Matrix::triu((W + Matrix::t(W)) / 2, k = 1L)
    el <- Matrix::summary(Ws)
    net <- data.frame(from = ids[el$i], to = ids[el$j], weight = el$x)

    st <- analyzeData(X, edgeStatsParam(weighted = TRUE), network = net)
    res <- if (identical(method, "moran")) .moran_test(st) else .geary_test(st)
    res
}

# binSpect(implementation = "fast"): the binSpectParam composition, through the
# giotto method. Takes binSpect's own arguments; the kmeans tuning arguments do
# not apply (the 2-means split is exact).
.binspect_fast <- function(gobject, spat_unit, feat_type, bin_method,
    expression_values, subset_feats, spatial_network_name,
    spatial_network_k, percentage_rank, adjust_method, return_gobject) {
    if (!is.null(spatial_network_k)) {
        stop("[binSpect] implementation = \"fast\" runs on one network; ",
            "build the kNN network and pass its spatial_network_name.",
            call. = FALSE)
    }
    if (!is.null(subset_feats)) {
        stop("[binSpect] implementation = \"fast\" tests every feature; ",
            "subset the result instead of `subset_feats`.", call. = FALSE)
    }
    param <- spatAutoCorParam("binspect",
        bin_method = bin_method[1L],
        percentage_rank = percentage_rank,
        adjust_method = adjust_method)
    analyzeData(gobject, param,
        spat_unit = spat_unit, feat_type = feat_type,
        expression_values = expression_values[1L],
        spatial_network_name = spatial_network_name,
        return_gobject = return_gobject)
}
