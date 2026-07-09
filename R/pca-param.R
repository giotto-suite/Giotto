#' @include dimension_reduction.R
NULL

# ============================================================================
# PCA parameter classes — extend GiottoClass::reduceParam so the
# decomposition is dispatched via reduceData(x, param). Distinct from
# processData (transforms data, same shape) and analyzeData (summary
# stats): reduceData returns a list(u, d, v, sdev, eigenvalues). Default
# methods call the existing .run_pca_biocsingular helper. Streaming
# backends (parquetExprStore in GiottoDisk) provide their own setMethod
# for randomPcaParam — Halko-style randomized SVD with streaming
# Cholesky-QR.
# ============================================================================

# ---- VIRTUAL + concrete ----------------------------------------------------

#' @rdname process_param
#' @exportClass pcaParam
setClass("pcaParam", contains = c("VIRTUAL", "reduceParam"))

#' @rdname process_param
#' @exportClass irlbaPcaParam
setClass("irlbaPcaParam", contains = "pcaParam")

#' @rdname process_param
#' @exportClass exactPcaParam
setClass("exactPcaParam", contains = "pcaParam")

#' @rdname process_param
#' @exportClass randomPcaParam
setClass("randomPcaParam", contains = "pcaParam")

# Auto sentinel: routes to a substrate-recommended concrete flavor at
# reduceData dispatch. Users construct via `pcaParam("auto", ...)`; the
# `reduceData(x, autoPcaParam)` method resolves to a concrete flavor
# based on `class(x)` and either dispatches (default) or returns the
# resolved param (`dry_run = TRUE`).
#' @rdname process_param
#' @exportClass autoPcaParam
setClass("autoPcaParam", contains = "pcaParam")


# ---- Factory ---------------------------------------------------------------

#' @rdname process_param
#' @title PCA parameter factory
#' @description
#' Construct a `pcaParam` for use with [reduceData()]. `runPCA()` builds
#' these internally; direct use is only needed when computing PCA on a
#' standalone matrix.
#'
#' Returns a list with `u`, `d`, `v`, and `sdev` when passed to
#' `reduceData()`.
#'
#' @param method one of `"auto"`, `"random"`, `"irlba"`, `"exact"`.
#'   `"auto"` defers the choice to the substrate: in-memory backends
#'   resolve to `"irlba"`; streaming backends (`parquetExprStore` in
#'   GiottoDisk) resolve to `"random"` / `"gram-eigen"` per their own
#'   heuristics.
#' @param ncp number of components. Default `50`.
#' @param center logical. Center columns. Default `TRUE`.
#' @param scale logical. Scale columns by SD. Default `TRUE`. Streaming
#'   backends (parquetExprStore) ignore this — scaling densifies the
#'   matrix and is incompatible with O(N*k) streaming.
#' @param feats_to_use character vector of feature IDs to subset to before
#'   PCA. `NULL` means all features. Useful for HVG selection.
#' @param n_oversamples integer. Halko oversampling parameter for `random`.
#'   Default `10`.
#' @param n_power_iter integer. Halko power iterations for `random`.
#'   Default `2`.
#' @param set_seed logical. Default `TRUE`.
#' @param seed_number integer. Default `1234`.
#' @param dry_run logical. Only meaningful when `method = "auto"`. When
#'   `TRUE`, `reduceData(x, param)` returns the substrate-resolved
#'   concrete `pcaParam` instead of running PCA. Useful for inspecting
#'   or testing the substrate's method selection. Default `FALSE`.
#' @param ... reserved.
#' @return Concrete `pcaParam` subclass object.
#' @examples
#' p <- pcaParam("random", ncp = 30)
#' p_auto <- pcaParam("auto", ncp = 30)  # substrate picks the flavor
#' # Inspect the substrate's choice without running PCA:
#' # reduceData(x, pcaParam("auto", dry_run = TRUE))
#' @export
pcaParam <- function(method        = c("auto", "random", "irlba", "exact"),
                      ncp           = 50L,
                      center        = TRUE,
                      scale         = TRUE,
                      feats_to_use  = NULL,
                      n_oversamples = 10L,
                      n_power_iter  = 2L,
                      set_seed      = TRUE,
                      seed_number   = 1234L,
                      dry_run       = FALSE,
                      ...) {
    method <- match.arg(tolower(method),
        c("auto", "random", "irlba", "exact"))
    cls <- switch(method,
        "auto"   = "autoPcaParam",
        "random" = "randomPcaParam",
        "irlba"  = "irlbaPcaParam",
        "exact"  = "exactPcaParam")
    p <- new(cls, param = list(...))
    p$method        <- method
    p$ncp           <- as.integer(ncp)
    p$center        <- isTRUE(center)
    p$scale         <- isTRUE(scale)
    p$feats_to_use  <- feats_to_use
    p$n_oversamples <- as.integer(n_oversamples)
    p$n_power_iter  <- as.integer(n_power_iter)
    p$set_seed      <- isTRUE(set_seed)
    p$seed_number   <- as.integer(seed_number)
    # dry_run is autoPcaParam-only: concrete flavors don't have this
    # concept (they always execute). Keeping it out of concrete flavors
    # avoids ambiguity when a resolved param is inspected downstream.
    if (identical(method, "auto")) p$dry_run <- isTRUE(dry_run)
    p
}


# ---- reduceData for autoPcaParam -----------------------------------------
# Each substrate registers its own method here to encode the routing
# choice. `allMatrix` covers the in-memory family (matrix / Matrix /
# DelayedArray / IterableMatrix) with IRLBA. Substrates in other
# packages (GiottoDisk::parquetExprStore) register their own method.
# The `ANY` fallback catches unregistered substrates with a warning +
# IRLBA -- last-resort routing that surfaces the "missing method"
# rather than a cryptic dispatch error, while still trying something
# sensible. Set `dry_run = TRUE` on the param to return the resolved
# concrete pcaParam instead of running PCA.

#' @rdname reduceData
setMethod("reduceData",
    signature(x = "allMatrix", param = "autoPcaParam"),
    function(x, param, ...) {
        resolved <- pcaParam(
            method        = "irlba",
            ncp           = param$ncp,
            center        = param$center,
            scale         = param$scale,
            feats_to_use  = param$feats_to_use,
            n_oversamples = param$n_oversamples,
            n_power_iter  = param$n_power_iter,
            set_seed      = param$set_seed,
            seed_number   = param$seed_number
        )
        if (isTRUE(param$dry_run)) return(resolved)
        reduceData(x, resolved, ...)
    }
)

#' @rdname reduceData
setMethod("reduceData",
    signature(x = "ANY", param = "autoPcaParam"),
    function(x, param, ...) {
        warning("[reduceData(", class(x)[1L], ", autoPcaParam)] no ",
            "substrate-specific auto routing registered for '",
            class(x)[1L], "'; falling back to method = \"irlba\". ",
            "Register a `reduceData(<class>, autoPcaParam)` method to ",
            "silence this warning.", call. = FALSE)
        resolved <- pcaParam(
            method        = "irlba",
            ncp           = param$ncp,
            center        = param$center,
            scale         = param$scale,
            feats_to_use  = param$feats_to_use,
            n_oversamples = param$n_oversamples,
            n_power_iter  = param$n_power_iter,
            set_seed      = param$set_seed,
            seed_number   = param$seed_number
        )
        if (isTRUE(param$dry_run)) return(resolved)
        reduceData(x, resolved, ...)
    }
)


# ---- Default methods on allMatrix ------------------------------------------
# Wraps the existing .run_pca_biocsingular helper. Returns a list with the
# fields runPCA already extracts from BiocSingular's output, plus sdev for
# streaming-backend symmetry.

.pca_default_method <- function(x, param, BSPARAM) {
    # Subset features if requested
    if (!is.null(param$feats_to_use)) {
        keep <- intersect(rownames(x), param$feats_to_use)
        if (length(keep) == 0L) {
            stop("[reduceData(allMatrix, pcaParam)] feats_to_use ",
                 "had zero overlap with rownames(x).", call. = FALSE)
        }
        x <- x[keep, , drop = FALSE]
    }
    # Input convention: x is gene x cell (Bioconductor / Giotto).
    # .run_pca_biocsingular's dense path expects cell x gene -> transpose.
    # Its IterableMatrix path (BPCells::svds) expects gene x cell -> no
    # transpose. Pick correctly so eigenvalues = d^2/(n_cells-1) for both.
    x_for_pca <- if (inherits(x, "IterableMatrix")) x else t_flex(x)
    res <- .run_pca_biocsingular(
        x = x_for_pca,
        ncp = param$ncp,
        center = param$center,
        scale = param$scale,
        BSPARAM = BSPARAM,
        set_seed = param$set_seed,
        seed_number = param$seed_number
    )
    # .run_pca_biocsingular returns eigenvalues = sdev^2 = d^2/(n-1)
    # where n is the number of observations (cells). Recover canonical
    # singular values d so streaming-backend output is directly comparable.
    n_obs <- ncol(x)   # x is gene x cell, so ncol = n_cells (observations)
    sdev <- sqrt(res$eigenvalues)
    d    <- sdev * sqrt(max(n_obs - 1L, 1L))
    list(
        u           = res$coords,    # cells x k embeddings (u * d)
        d           = d,              # singular values (matches streaming)
        v           = res$loadings,  # genes x k loadings
        sdev        = sdev,
        eigenvalues = res$eigenvalues
    )
}

#' @rdname reduceData
setMethod("reduceData",
    signature(x = "allMatrix", param = "irlbaPcaParam"),
    function(x, param, ...) .pca_default_method(x, param, BSPARAM = "irlba")
)

#' @rdname reduceData
setMethod("reduceData",
    signature(x = "allMatrix", param = "exactPcaParam"),
    function(x, param, ...) .pca_default_method(x, param, BSPARAM = "exact")
)

#' @rdname reduceData
setMethod("reduceData",
    signature(x = "allMatrix", param = "randomPcaParam"),
    function(x, param, ...) .pca_default_method(x, param, BSPARAM = "random")
)


# ---- BiocSingularParam pass-through ---------------------------------------
# Users already comfortable with BiocSingular's `runPCA(BSPARAM = ...)` API
# can pass a pre-built `IrlbaParam()` / `ExactParam()` / `RandomParam()` /
# `FastAutoParam()` directly, exposing knobs (`deferred`, `fold`, ...) that
# Giotto's pcaParam wrapper doesn't surface. Bioc-native dispatch only --
# IterableMatrix and other non-dense Giotto-native backends still need
# their own setMethod or the pcaParam factory.

#' @rdname reduceData
setMethod("reduceData",
    signature(x = "allMatrix", param = "BiocSingularParam"),
    function(x, param,
             ncp           = 50L,
             center        = TRUE,
             scale         = TRUE,
             feats_to_use  = NULL,
             set_seed      = TRUE,
             seed_number   = 1234L,
             ...) {
        if (!is.null(feats_to_use)) {
            keep <- intersect(rownames(x), feats_to_use)
            if (length(keep) == 0L) {
                stop("[reduceData(allMatrix, BiocSingularParam)] feats_to_use ",
                     "had zero overlap with rownames(x).", call. = FALSE)
            }
            x <- x[keep, , drop = FALSE]
        }
        x_for_pca <- if (inherits(x, "IterableMatrix")) x else t_flex(x)
        res <- .run_pca_biocsingular(
            x           = x_for_pca,
            ncp         = ncp,
            center      = center,
            scale       = scale,
            BSPARAM     = param,
            set_seed    = set_seed,
            seed_number = seed_number
        )
        n_obs <- ncol(x)
        sdev <- sqrt(res$eigenvalues)
        d    <- sdev * sqrt(max(n_obs - 1L, 1L))
        list(
            u           = res$coords,
            d           = d,
            v           = res$loadings,
            sdev        = sdev,
            eigenvalues = res$eigenvalues
        )
    }
)
