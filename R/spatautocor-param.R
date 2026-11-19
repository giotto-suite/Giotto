#' @include spatautocor-stats.R
NULL

# ============================================================================
# Spatial autocorrelation parameter classes. A VIRTUAL spatAutoCorParam with
# one concrete class per statistic, built by spatAutoCorParam(). Each runs as
# a composition of building blocks -- thresholdParam("binarize") (binSpect
# only), then GiottoClass's edgeStatsParam -- followed by the statistic in
# spatautocor-stats.R, so a backend that implements those two blocks gets all
# three statistics without code of its own.
# ============================================================================

#' @title Spatial autocorrelation parameter classes
#' @name spatAutoCorParam-class
#' @aliases moranParam-class gearyParam-class binSpectParam-class
#' @description
#' Parameter classes for spatial autocorrelation of features over a spatial
#' network: `moranParam` (Moran's I), `gearyParam` (Geary's C) and
#' `binSpectParam` (binSpect, a join count on binarized features). All extend
#' the virtual `spatAutoCorParam`. Build them with [spatAutoCorParam()].
#' @exportClass spatAutoCorParam
setClass("spatAutoCorParam", contains = c("VIRTUAL", "analyzeParam"))

#' @rdname spatAutoCorParam-class
#' @exportClass moranParam
setClass("moranParam", contains = "spatAutoCorParam")

#' @rdname spatAutoCorParam-class
#' @exportClass gearyParam
setClass("gearyParam", contains = "spatAutoCorParam")

#' @rdname spatAutoCorParam-class
#' @exportClass binSpectParam
setClass("binSpectParam", contains = "spatAutoCorParam")


#' @title Spatial autocorrelation parameters
#' @name spatAutoCorParam
#' @description
#' Settings for detecting spatially autocorrelated features over a spatial
#' network. Significance comes from the randomisation null in closed form:
#' the exact mean and variance of the statistic over all permutations of a
#' feature's values across cells (Cliff and Ord 1981), so nothing is shuffled.
#'
#' * `"moran"` — Moran's I of each feature's values.
#' * `"geary"` — Geary's C of each feature's values.
#' * `"binspect"` — binSpect: each feature is binarized per cell (see
#'   [threshold_binarize]) and tested with the black-black join
#'   count, the number of edges with both ends on. This is Moran's I of the
#'   binarized feature up to the network's degree variation, and the 2x2
#'   table of edge types it summarises is the one binSpect has always used
#'   (Dries et al. 2021); the odds ratio of that table is reported alongside.
#' @param method one of `"moran"`, `"geary"`, `"binspect"`; `"joincount"` is
#'   another name for `"binspect"` (the same test, the same result columns).
#' @param null `"permutation"` (default) or `"covariate"`. The covariate null
#'   asks whether a feature clusters beyond what a cell grouping explains:
#'   cells are grouped by `covariate`, and Moran's I is taken of each value
#'   minus its group's mean, or binSpect's join count is compared with cells
#'   switching on at their group's rate. Use it when something other than
#'   location -- sequencing depth, cell type -- varies across the tissue and
#'   would otherwise make many features look spatial. Not available for
#'   `"geary"`; binSpect needs an unweighted network.
#' @param covariate for `null = "covariate"` on a `giotto` object: `"libsize"`
#'   (the cells' total raw counts) or a cell metadata column. A numeric
#'   covariate is cut into `n_groups` quantile groups; a categorical one is
#'   used as is. Matrix and store inputs take the grouping directly as
#'   `groups =`.
#' @param n_groups number of quantile groups for a numeric covariate.
#' @param weighted logical. Use the network's edge weights.
#' @param adjust_method p-value adjustment, passed to [stats::p.adjust()].
#' @param bin_method,percentage_rank binarization for `"binspect"`: the
#'   `threshold` rule and `percentage_rank` of [threshold_binarize].
#' @param ... additional named entries to attach to `@param`.
#' @returns A concrete `spatAutoCorParam` subclass object.
#' @section Features that do not vary:
#' A feature with the same value in every cell (all zero, or constant), or
#' for binSpect one binarized all on or all off, has no randomisation null.
#' Its `z`, `p.value` and `adj.p.value` are `NA`, as is a statistic that is
#' undefined for it (Moran's I, Geary's C); it is left out of the p-value
#' adjustment and sorts last.
#' @references
#' Cliff, A. D. and Ord, J. K. (1981). Spatial Processes: Models and
#' Applications. Pion.
#'
#' Dries, R. et al. (2021). Giotto: a toolbox for integrative analysis and
#' visualization of spatial expression data. Genome Biology 22, 78.
#' @examples
#' p <- spatAutoCorParam("moran")
#' p_cov <- spatAutoCorParam("moran", null = "covariate", covariate = "libsize")
#' p_bin <- spatAutoCorParam("binspect", bin_method = "rank")
#' @export
spatAutoCorParam <- function(method = c("moran", "geary", "binspect", "joincount"),
    null = c("permutation", "covariate"),
    covariate = "libsize",
    n_groups = 10L,
    weighted = FALSE,
    adjust_method = "fdr",
    bin_method = c("kmeans", "rank"),
    percentage_rank = 30,
    ...) {
    method <- match.arg(method, c("moran", "geary", "binspect", "joincount"))
    if (identical(method, "joincount")) method <- "binspect"
    null <- match.arg(null, c("permutation", "covariate"))
    if (null == "covariate" && method == "geary") {
        stop("[spatAutoCorParam] the covariate null is available for ",
            "\"moran\" and \"binspect\".", call. = FALSE)
    }
    checkmate::assert_string(covariate)
    checkmate::assert_count(n_groups, positive = TRUE)
    checkmate::assert_flag(weighted)
    adjust_method <- match.arg(adjust_method, stats::p.adjust.methods)
    cls <- switch(method,
        moran = "moranParam",
        geary = "gearyParam",
        binspect = "binSpectParam"
    )
    p <- new(cls, param = list(...))
    p$method <- method
    p$null <- null
    if (null == "covariate") {
        p$covariate <- covariate
        p$n_groups <- as.integer(n_groups)
    }
    p$weighted <- weighted
    p$adjust_method <- adjust_method
    if (identical(method, "binspect")) {
        p$bin_method <- match.arg(bin_method, c("kmeans", "rank"))
        p$percentage_rank <- percentage_rank
    }
    p
}
