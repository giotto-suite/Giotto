# deconvParam ####
#
# Cell type deconvolution as an `analyzeData` verb, the sibling of
# `enrichParam` in enrich-param.R: a VIRTUAL parent under `analyzeParam`, one
# concrete class per method, a factory, and an `analyzeData` method that
# dispatches on the *expression object*. A disk-backed store therefore arrives
# intact and can route to its own streaming method (GiottoDisk) with no change
# here.
#
# The algorithm is not in this package: it is the spatialDWLS package
# (github: drieslab/spatialDWLS), which Giotto suggests, the way the motif
# layer uses smotif. Giotto owns the param class, the verb and
# the gobject plumbing in runDWLSDeconv(); the `analyzeData` method below is
# the only code here that calls the engine.

#' @title Deconvolution params
#' @name deconv_param
#' @description
#' Parameter objects for spatial cell type deconvolution. Pass one to
#' [analyzeData()] together with a signature matrix and a cluster vector, or
#' use the [runDWLSDeconv()] wrapper, which builds the param for you.
#'
#' @param method character. Only `"DWLS"` for now.
#' @param ... method-specific parameters; see Details.
#' @details
#' DWLS:
#' \describe{
#'   \item{`logbase`}{numeric. log base of the expression values (default 2).}
#'   \item{`n_cell`}{numeric. expected number of cells per spot (default 50).
#'     A cell type is kept for a spot in the second round only when its first
#'     round proportion is at least `1 / n_cell`.}
#'   \item{`cutoff`}{numeric. PAGE z-score a cell type must reach somewhere in
#'     a cluster to be a candidate for that cluster (default 2).}
#'   \item{`na_as_zero`}{logical. spots that cannot be solved (no signal on
#'     the signature features) are `NA` by default; `TRUE` returns them as 0.}
#'   \item{`BPPARAM`}{a `BiocParallelParam` for the per-spot solves. `NULL`
#'     (default) uses `getOption("giotto.bpparam")` when that is a parallel
#'     param, and runs serially otherwise.}
#'   \item{`chunk_size`}{integer. spots per parallel task (default 256).}
#' }
#' @md
#' @returns a [deconvParam-class]-inheriting object
#' @examples
#' p <- deconvParam("DWLS", n_cell = 20)
#' p$n_cell
#' @family spatial deconvolution
#' @seealso [analyze_param], [runDWLSDeconv()]
NULL

#' @rdname deconv_param
#' @exportClass deconvParam
setClass("deconvParam", contains = c("VIRTUAL", "analyzeParam"))

#' @rdname deconv_param
#' @exportClass dwlsDeconvParam
setClass("dwlsDeconvParam", contains = "deconvParam")


# param factory ####

#' @rdname deconv_param
#' @export
deconvParam <- function(method = "DWLS", ...) {
    method <- match.arg(toupper(method), "DWLS")
    switch(method,
        "DWLS" = .deconv_param_dwls(...)
    )
}

#' @keywords internal
#' @noRd
.deconv_param_dwls <- function(...) {
    p <- new("dwlsDeconvParam", param = list(...))
    p$logbase <- as.numeric(p$logbase %null% 2)
    p$n_cell <- as.numeric(p$n_cell %null% 50)
    p$cutoff <- as.numeric(p$cutoff %null% 2)
    p$na_as_zero <- isTRUE(p$na_as_zero)
    p$chunk_size <- as.integer(p$chunk_size %null% 256L)
    if (!is.null(p$BPPARAM)) {
        checkmate::assert_class(p$BPPARAM, "BiocParallelParam")
    }
    if (length(p$n_cell) != 1L || !is.finite(p$n_cell) || p$n_cell <= 0) {
        stop("[deconvParam] `n_cell` must be a single positive number.",
            call. = FALSE)
    }
    p
}


# DWLS ####

#' @title spatialDWLS deconvolution
#' @name deconv_dwls
#' @description
#' The [analyzeData()] method for a [dwlsDeconvParam-class]. The computation
#' is done by the \pkg{spatialDWLS} package (`spatialDWLS::dwls()`), which is
#' suggested rather than imported: install it with
#' `remotes::install_github("drieslab/spatialDWLS")`.
#' @param x expression values, log-scaled, features by spots. A `matrix`, a
#'   `Matrix`, or anything an engine registered against
#'   [dwlsDeconvParam-class] accepts.
#' @param param a [dwlsDeconvParam-class].
#' @inheritParams deconvolution_params
#' @param cell_clusters cluster assignment of each spot. Either named by
#'   spot ID (matched against `colnames(x)`) or unnamed and in the same order
#'   as the columns of `x`.
#' @param ... unused. Declared because [analyzeData()] carries it, so an
#'   engine registered from another package can accept its own arguments.
#' @md
#' @returns a `data.table` of `cell_ID` and one proportion column per cell
#'   type. Each solved spot sums to one; unsolvable spots are `NA` (or 0, see
#'   `na_as_zero` in [deconvParam()]).
#' @family spatial deconvolution
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' expr <- GiottoClass::getExpression(g, values = "normalized",
#'                                    output = "matrix")
#' md <- GiottoClass::pDataDT(g)
#' sm <- makeSignMatrixDWLS(g,
#'     sign_gene = rownames(expr)[1:100],
#'     cell_type_vector = as.character(md$leiden_clus)
#' )
#' clus <- stats::setNames(md$leiden_clus, md$cell_ID)
#'
#' res <- analyzeData(expr, deconvParam("DWLS"),
#'     sign_matrix = sm, cell_clusters = clus)
#' head(res)
#' @export
setMethod("analyzeData",
    signature(x = "ANY", param = "dwlsDeconvParam"),
    function(x, param, ..., sign_matrix, cell_clusters) {
        .dwls_engine_check()
        prop <- spatialDWLS::dwls(
            expr = x,
            sign_matrix = sign_matrix,
            clusters = cell_clusters,
            n_cell = param$n_cell,
            cutoff = param$cutoff,
            logbase = param$logbase,
            na_as_zero = param$na_as_zero,
            BPPARAM = .dwls_bpparam(param$BPPARAM),
            chunk_size = param$chunk_size
        )
        out <- data.table::data.table(cell_ID = colnames(prop))
        .enrich_check_contract(
            cbind(out, data.table::as.data.table(t(prop))),
            "analyzeData(dwlsDeconvParam)"
        )
    }
)

# The engine is suggested, not imported: stop with the install line when it
# is missing.
#' @keywords internal
#' @noRd
.dwls_engine_check <- function() {
    package_check(
        pkg_name = "spatialDWLS",
        repository = "github:drieslab/spatialDWLS"
    )
}

# The param's own BPPARAM, else the session's `giotto.bpparam`. The engine runs
# serially for NULL or a SerialParam, so neither warns the way lapply_flex()
# does on every serial run.
#' @keywords internal
#' @noRd
.dwls_bpparam <- function(BPPARAM = NULL) {
    BPPARAM %null% getOption("giotto.bpparam")
}
