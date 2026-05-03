#' @include variable_genes.R
NULL

# ============================================================================
# HVG parameter classes — extend processParam family. Same dispatch pattern
# as norm/scale/qc/filter. The mathematical helpers (.calc_expr_cov_stats,
# .calc_cov_group_hvf, .calc_cov_loess_hvf, .calc_var_hvf) are reused
# verbatim in the default methods so dgCMatrix / IterableMatrix backends
# produce bit-for-bit identical numerical output.
# ============================================================================

# ---- VIRTUAL + concrete ----------------------------------------------------

#' @rdname process_param
#' @exportClass hvgParam
setClass("hvgParam", contains = c("VIRTUAL", "processParam"))

#' @rdname process_param
#' @exportClass covGroupsHvgParam
setClass("covGroupsHvgParam", contains = "hvgParam")

#' @rdname process_param
#' @exportClass covLoessHvgParam
setClass("covLoessHvgParam", contains = "hvgParam")

#' @rdname process_param
#' @exportClass varPResidHvgParam
setClass("varPResidHvgParam", contains = "hvgParam")


# ---- Factory ---------------------------------------------------------------

#' @rdname process_param
#' @title HVG parameter factory
#' @description
#' Construct an `hvgParam` for use with [processData()]. `calculateHVF()`
#' builds these internally; direct use is only needed when computing HVG
#' on a standalone matrix.
#'
#' Returns a `data.table` of feature-level statistics with a `selected`
#' column ("yes"/"no") when passed to `processData()`.
#'
#' @param method one of `"cov_groups"`, `"cov_loess"`, `"var_p_resid"`.
#' @param expression_threshold detection threshold passed to
#'   `.calc_expr_cov_stats`. Default `0`.
#' @param nr_expression_groups (`cov_groups`) number of mean-expression
#'   bins. Default `20`.
#' @param zscore_threshold (`cov_groups`) within-bin z-score cutoff.
#'   Default `1.5`.
#' @param difference_in_cov (`cov_loess`) cov - predicted_cov cutoff.
#'   Default `0.1`.
#' @param var_threshold (`var_p_resid`) variance cutoff. Default `1.5`.
#' @param var_number (`var_p_resid`) keep top-N by variance. Default
#'   `NULL` (use threshold instead).
#' @param calc_gini logical. Whether to compute Gini index. Default
#'   `FALSE` (faster, fine for HVG selection).
#' @param ... reserved.
#' @return Concrete `hvgParam` subclass object.
#' @examples
#' p <- hvgParam("cov_loess", difference_in_cov = 0.1)
#' @export
hvgParam <- function(method = c("cov_groups", "cov_loess", "var_p_resid"),
                      expression_threshold = 0,
                      nr_expression_groups = 20,
                      zscore_threshold     = 1.5,
                      difference_in_cov    = 0.1,
                      var_threshold        = 1.5,
                      var_number           = NULL,
                      calc_gini            = FALSE,
                      ...) {
    method <- match.arg(tolower(method),
        c("cov_groups", "cov_loess", "var_p_resid"))
    cls <- switch(method,
        "cov_groups"  = "covGroupsHvgParam",
        "cov_loess"   = "covLoessHvgParam",
        "var_p_resid" = "varPResidHvgParam")
    p <- new(cls, param = list(...))
    p$expression_threshold <- p$expression_threshold %null% expression_threshold
    p$nr_expression_groups <- p$nr_expression_groups %null% nr_expression_groups
    p$zscore_threshold     <- p$zscore_threshold     %null% zscore_threshold
    p$difference_in_cov    <- p$difference_in_cov    %null% difference_in_cov
    p$var_threshold        <- p$var_threshold        %null% var_threshold
    p$var_number           <- p$var_number           %null% var_number
    p$calc_gini            <- p$calc_gini            %null% calc_gini
    p
}


# ---- Default methods on allMatrix ------------------------------------------
# Each lifts the body that was inline in calculateHVF's switch() statement,
# without the plot generation (orchestrated at calculateHVF level).

#' @rdname processData
setMethod("processData",
    signature(x = "allMatrix", param = "covGroupsHvgParam"),
    function(x, param, ...) {
        stats <- .calc_expr_cov_stats(
            x,
            expression_threshold = param$expression_threshold,
            calc_gini            = param$calc_gini
        )
        res <- .calc_cov_group_hvf(
            stats,
            nr_expression_groups = param$nr_expression_groups,
            zscore_threshold     = param$zscore_threshold,
            show_plot   = FALSE, return_plot = FALSE, save_plot = FALSE
        )
        res$dt
    }
)

#' @rdname processData
setMethod("processData",
    signature(x = "allMatrix", param = "covLoessHvgParam"),
    function(x, param, ...) {
        stats <- .calc_expr_cov_stats(
            x,
            expression_threshold = param$expression_threshold,
            calc_gini            = param$calc_gini
        )
        res <- .calc_cov_loess_hvf(
            stats,
            difference_in_cov = param$difference_in_cov,
            show_plot   = FALSE, return_plot = FALSE, save_plot = FALSE
        )
        res$dt
    }
)

#' @rdname processData
setMethod("processData",
    signature(x = "allMatrix", param = "varPResidHvgParam"),
    function(x, param, ...) {
        res <- .calc_var_hvf(
            scaled_matrix = x,
            var_threshold = param$var_threshold,
            var_number    = param$var_number,
            show_plot   = FALSE, return_plot = FALSE, save_plot = FALSE,
            use_parallel = FALSE
        )
        res$dt
    }
)
