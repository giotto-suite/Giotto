#' @include filter.R
NULL

# ============================================================================
# Filter parameter class — extends GiottoClass::filterParam so the cell +
# feature mask computation goes through the filterData(x, param) S4
# dispatch. Distinct from processData (transforms data) and analyzeData
# (computes summary stats): filterData returns a selection
# (list(feats_keep, cells_keep) of character IDs). The mask APPLICATION
# step (subsetGiotto) stays in filterGiotto and is unchanged for in-memory
# backends.
# ============================================================================

# ---- Concrete --------------------------------------------------------------

#' @rdname process_param
#' @exportClass defaultFilterParam
setClass("defaultFilterParam", contains = "filterParam")


# ---- Factory ---------------------------------------------------------------

#' @rdname process_param
#' @title Filter parameter factory
#' @description
#' Construct a `filterParam` carrying the cell + feature filter thresholds
#' used by [filterGiotto()]. `filterGiotto` builds this internally; direct
#' use is only needed when computing masks on a standalone matrix.
#'
#' Returns a `list(feats_keep, cells_keep)` of character ID vectors when
#' passed to [filterData()].
#'
#' @param method only `"default"` for now (Giotto's two-stage filter).
#' @param expression_threshold numeric. A value `>= expression_threshold`
#'   counts as detected. Default `1`.
#' @param feat_det_in_min_cells integer. Keep features detected in at least
#'   this many cells. Default `100`.
#' @param min_det_feats_per_cell integer. Keep cells with at least this many
#'   detected features (counted only over kept features — Giotto's
#'   two-stage convention). Default `100`.
#' @param ... reserved.
#' @return A `defaultFilterParam` object.
#' @examples
#' p <- filterParam(expression_threshold = 1,
#'                   feat_det_in_min_cells = 5,
#'                   min_det_feats_per_cell = 3)
#' @export
filterParam <- function(method = "default",
                         expression_threshold   = 1,
                         feat_det_in_min_cells  = 100,
                         min_det_feats_per_cell = 100,
                         ...) {
    method <- match.arg(tolower(method), c("default"))
    cls <- switch(method, "default" = "defaultFilterParam")
    p <- new(cls, param = list(...))
    p$expression_threshold   <- p$expression_threshold   %null% expression_threshold
    p$feat_det_in_min_cells  <- p$feat_det_in_min_cells  %null% feat_det_in_min_cells
    p$min_det_feats_per_cell <- p$min_det_feats_per_cell %null% min_det_feats_per_cell
    p
}


# ---- Default method on allMatrix -------------------------------------------
# Body is the exact two-stage logic that lives inline in filterGiotto today,
# so dgCMatrix / Matrix / IterableMatrix backends produce bit-for-bit
# identical mask output before vs after the refactor.

#' @rdname filterData
setMethod("filterData",
    signature(x = "allMatrix", param = "filterParam"),
    function(x, param, ...) {
        thr   <- param$expression_threshold
        f_min <- param$feat_det_in_min_cells
        c_min <- param$min_det_feats_per_cell

        # 1. feature mask
        filter_index_feats <- rowSums_flex(x >= thr) >= f_min
        selected_feat_ids <- names(filter_index_feats[filter_index_feats == TRUE])

        # 2. cell mask (only on kept features)
        filter_index_cells <- colSums_flex(x[filter_index_feats, ] >= thr) >= c_min
        selected_cell_ids <- names(filter_index_cells[filter_index_cells == TRUE])

        list(feats_keep = selected_feat_ids,
             cells_keep = selected_cell_ids)
    }
)

#' @rdname filterData
setMethod("filterData",
    signature(x = "IterableMatrix", param = "filterParam"),
    function(x, param, ...) {
        thr   <- param$expression_threshold
        f_min <- param$feat_det_in_min_cells
        c_min <- param$min_det_feats_per_cell

        filter_index_feats <- rowSums_flex(x >= thr) >= f_min
        selected_feat_ids <- names(filter_index_feats[filter_index_feats == TRUE])

        filter_index_cells <- colSums_flex(x[filter_index_feats, ] >= thr) >= c_min
        selected_cell_ids <- names(filter_index_cells[filter_index_cells == TRUE])

        list(feats_keep = selected_feat_ids,
             cells_keep = selected_cell_ids)
    }
)
