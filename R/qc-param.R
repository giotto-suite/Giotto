#' @include normalize.R
NULL

# ============================================================================
# QC parameter classes — extend the processParam family so QC computations
# go through the same processData(x, param) S4 dispatch used by normalize /
# scale / adjust / threshold. This lets backends (parquetExprStore,
# IterableMatrix, dgCMatrix, ...) provide their own implementations
# without changing addCellStatistics / addFeatStatistics signatures.
# ============================================================================

# ---- VIRTUAL classes -------------------------------------------------------

#' @rdname process_param
#' @exportClass qcParam
setClass("qcParam", contains = c("VIRTUAL", "processParam"))

#' @rdname process_param
#' @exportClass cellQcParam
setClass("cellQcParam", contains = "qcParam")

#' @rdname process_param
#' @exportClass featQcParam
setClass("featQcParam", contains = "qcParam")


# ---- Factory ---------------------------------------------------------------

#' @rdname process_param
#' @title QC parameter factories
#' @description
#' Construct a `qcParam` (cell- or feature-level) for use with
#' [processData()]. `addCellStatistics()` / `addFeatStatistics()` /
#' `addStatistics()` build these internally when called on a giotto
#' object, so direct use is only needed when computing QC stats on a
#' standalone matrix.
#'
#' @param level `"cell"` or `"feat"`.
#' @param detection_threshold numeric. A feature is detected in a cell if
#'   the value is strictly greater than this. Default `0`.
#' @param ... reserved for future params.
#' @return A `cellQcParam` or `featQcParam` object.
#' @examples
#' p <- qcParam("cell", detection_threshold = 0)
#' @export
qcParam <- function(level = c("cell", "feat"),
                     detection_threshold = 0, ...) {
    level <- match.arg(tolower(level), c("cell", "feat"))
    cls <- switch(level, "cell" = "cellQcParam", "feat" = "featQcParam")
    p <- new(cls, param = list(...))
    p$detection_threshold <- p$detection_threshold %null% detection_threshold
    p
}


# ---- Default methods on allMatrix ------------------------------------------
# Bodies are exact lifts of the inline code currently in
# addCellStatistics / addFeatStatistics, so dgCMatrix / Matrix /
# IterableMatrix backends produce bit-for-bit identical numerical output
# pre- and post-refactor.

#' @rdname processData
setMethod("processData",
    signature(x = "allMatrix", param = "cellQcParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        data.table::data.table(
            cells      = colnames(x),
            nr_feats   = colSums_flex(x > thr),
            perc_feats = (colSums_flex(x > thr) / nrow(x)) * 100,
            total_expr = colSums_flex(x)
        )
    }
)

#' @rdname processData
setMethod("processData",
    signature(x = "allMatrix", param = "featQcParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        feat_stats <- data.table::data.table(
            feats      = rownames(x),
            nr_cells   = rowSums_flex(x > thr),
            perc_cells = (rowSums_flex(x > thr) / ncol(x)) * 100,
            total_expr = rowSums_flex(x),
            mean_expr  = rowMeans_flex(x)
        )
        mean_expr_det <- NULL  # NSE binding for data.table
        feat_stats[, mean_expr_det := .mean_expr_det_test(
            x, detection_threshold = thr
        )]
        feat_stats
    }
)

# Mirror methods on dgCMatrix / IterableMatrix delegate to allMatrix —
# allMatrix is the union (matrix, Matrix), so dgCMatrix is already covered.
# IterableMatrix needs an explicit method because it's NOT in the allMatrix
# union but rowSums_flex/colSums_flex/rowMeans_flex handle it.

#' @rdname processData
setMethod("processData",
    signature(x = "IterableMatrix", param = "cellQcParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        data.table::data.table(
            cells      = colnames(x),
            nr_feats   = colSums_flex(x > thr),
            perc_feats = (colSums_flex(x > thr) / nrow(x)) * 100,
            total_expr = colSums_flex(x)
        )
    }
)

#' @rdname processData
setMethod("processData",
    signature(x = "IterableMatrix", param = "featQcParam"),
    function(x, param, ...) {
        thr <- param$detection_threshold %null% 0
        feat_stats <- data.table::data.table(
            feats      = rownames(x),
            nr_cells   = rowSums_flex(x > thr),
            perc_cells = (rowSums_flex(x > thr) / ncol(x)) * 100,
            total_expr = rowSums_flex(x),
            mean_expr  = rowMeans_flex(x)
        )
        mean_expr_det <- NULL
        feat_stats[, mean_expr_det := .mean_expr_det_test(
            x, detection_threshold = thr
        )]
        feat_stats
    }
)
