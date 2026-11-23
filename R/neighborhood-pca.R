#' @include pca-param.R
NULL

# neighborhoodPCAParam ####
#
# A native implementation of the BANKSY embedding (Singhal et al. 2024) as a
# `reduceData()` verb (adr/0008): each cell's values are joined by the
# Gaussian-weighted mean of its k nearest spatial neighbours, both blocks are
# scaled, weighted by sqrt(1 - lambda) and sqrt(lambda), and reduced by PCA.
# The spatial kNN comes from GiottoClass `createNetwork()` and the PCA from
# `reduceData(x, pcaParam("irlba"))`; only the neighbour weights are new.
# GiottoConnect::runBanksy() runs the original Bioconductor package with the
# same verb and defaults, and is the reference.

#' @title Neighborhood PCA parameters
#' @name neighborhoodPCAParam
#' @description
#' Parameters for [reduceData()] that compute the BANKSY embedding of
#' Singhal et al. (2024): each cell's values are joined by the weighted mean
#' of its `k_geom` nearest spatial neighbours, the two blocks are scaled
#' per feature, weighted by \eqn{\sqrt{1 - \lambda}} and
#' \eqn{\sqrt{\lambda}}, and reduced by PCA. Neighbour weights are
#' \eqn{\exp(-d^2 / m^2)}, with \eqn{m} the median distance of the cell's
#' `k_geom` neighbours, normalized to sum to one per cell, as in the
#' \pkg{Banksy} package.
#'
#' With `compute_agf = TRUE` a third block adds BANKSY's azimuthal Gabor
#' filter (AGF): the magnitude of \eqn{\sum_i w_{ji} e^{i\phi_{ji}}
#' (x_i - \bar{x}_j)}, the gradient of the neighbourhood around cell
#' \eqn{j} (\eqn{\phi} the angle to neighbour \eqn{i}, \eqn{\bar{x}_j} the
#' mean of its neighbours), and \eqn{\lambda} is split over the blocks as
#' \eqn{\sqrt{1 - \lambda}}, \eqn{\sqrt{2\lambda/3}} and
#' \eqn{\sqrt{\lambda/3}}, as in \pkg{Banksy}.
#'
#' This is a reimplementation of BANKSY's embedding. **If you use it, please
#' cite the original paper** (see References).
#' `GiottoConnect::runBanksy()` runs the original \pkg{Banksy} package.
#' @param lambda numeric in `[0, 1]`, the weight of the neighbourhood. About
#'   0.2 refines cell types with their surroundings; about 0.8 finds spatial
#'   domains.
#' @param k_geom integer, number of spatial neighbours; with
#'   `compute_agf = TRUE` one value, or two (neighbours for the mean and for
#'   the gradient)
#' @param compute_agf logical, add the AGF (gradient) block
#' @param npcs integer, number of principal components
#' @param absolute logical. Take the absolute value of the neighbour means,
#'   as \pkg{Banksy} 1.6 does. This makes no difference for non-negative
#'   values such as expression, but discards the sign of signed inputs such
#'   as an embedding, so it is `FALSE` by default.
#' @param seed integer, seed of the PCA
#' @returns A `neighborhoodPCAParam` object.
#' @references Singhal V, Chou N, Lee J, et al. BANKSY unifies cell typing
#'   and tissue domain segmentation for scalable spatial omics data analysis.
#'   *Nature Genetics* 56, 431–441 (2024). \doi{10.1038/s41588-024-01664-3}
#' @seealso [runNeighborhoodPCA()], `GiottoConnect::runBanksy()`
#' @examples
#' p <- neighborhoodPCAParam(lambda = 0.8, k_geom = 15)
#' p$lambda
#' @md
#' @export
neighborhoodPCAParam <- function(lambda = 0.2, k_geom = 15L, npcs = 20L,
    compute_agf = FALSE, absolute = FALSE, seed = 1234L) {
    checkmate::assert_number(lambda, lower = 0, upper = 1)
    checkmate::assert_integerish(k_geom, lower = 1, any.missing = FALSE,
        min.len = 1L, max.len = 2L)
    checkmate::assert_flag(compute_agf)
    checkmate::assert_count(npcs, positive = TRUE)
    checkmate::assert_flag(absolute)
    checkmate::assert_count(seed)
    p <- new("neighborhoodPCAParam", param = list())
    p$lambda <- lambda
    p$k_geom <- as.integer(k_geom)
    p$npcs <- as.integer(npcs)
    p$compute_agf <- compute_agf
    p$absolute <- absolute
    p$seed <- as.integer(seed)
    p
}

#' @rdname neighborhoodPCAParam
#' @exportClass neighborhoodPCAParam
setClass("neighborhoodPCAParam", contains = "reduceParam")


# reduceData(allMatrix, neighborhoodPCAParam) ####

#' @rdname neighborhoodPCAParam
#' @param x features x cells values (`matrix` or `Matrix`)
#' @param param a `neighborhoodPCAParam`
#' @param ... unused
#' @param spatial_locs cells x 2 (or 3) coordinates, with the cell IDs as
#'   row names
#' @details `reduceData(x, param, spatial_locs = )` returns the list of
#'   [pcaParam()]'s `reduceData()` methods (`u`: cells x `npcs` coordinates,
#'   `d`, `v`: loadings of the joint features, `sdev`, `eigenvalues`).
#' @export
setMethod("reduceData",
    signature(x = "allMatrix", param = "neighborhoodPCAParam"),
    function(x, param, ..., spatial_locs) {
        if (is.null(colnames(x)) || is.null(rownames(spatial_locs))) {
            stop("[neighborhoodPCAParam] `x` needs cell IDs as column names ",
                "and `spatial_locs` as row names", call. = FALSE)
        }
        locs <- as.matrix(spatial_locs)[colnames(x), , drop = FALSE]
        if (anyNA(locs)) {
            stop("[neighborhoodPCAParam] some cells have no spatial location",
                call. = FALSE)
        }
        k0 <- param$k_geom[1L]
        x <- as.matrix(x)
        ops <- .neighbor_operators(locs, k0)
        nbr <- as.matrix(x %*% Matrix::t(ops$w))
        if (isTRUE(param$absolute)) nbr <- abs(nbr)
        rownames(nbr) <- paste0(rownames(x), "_H0")
        lam <- c(1 - param$lambda, param$lambda)
        blocks <- list(x, nbr)
        if (isTRUE(param$compute_agf)) {
            k1 <- param$k_geom[length(param$k_geom)]
            if (k1 != k0) ops <- .neighbor_operators(locs, k1)
            agf <- .neighbor_agf(x, ops)
            rownames(agf) <- paste0(rownames(x), "_H1")
            blocks <- c(blocks, list(agf))
            # Banksy's getLambdas(): harmonic m weighs lambda * 2^-m, rescaled
            lam <- c(1 - param$lambda, param$lambda * c(2, 1) / 3)
        }
        joint <- do.call(rbind, Map(function(l, b) sqrt(l) * .scale_rows(b),
            lam, blocks))
        joint <- joint[MatrixGenerics::rowSds(joint) != 0, , drop = FALSE]
        reduceData(joint, pcaParam("irlba", ncp = param$npcs, center = TRUE,
            scale = FALSE, seed_number = param$seed))
    }
)

# Neighbour operators of each cell's k nearest spatial neighbours (itself
# excluded), cells x cells, row j for cell j:
# w   weights exp(-d^2 / median(d)^2) normalized per cell (Banksy kNN_median)
# wc, ws  w * cos(phi), w * sin(phi), phi the angle from j to the neighbour
# a   1 / k: the plain neighbour mean (Banksy's fscale() centring)
.neighbor_operators <- function(locs, k) {
    distance <- from <- weight <- NULL # NSE
    ids <- rownames(locs)
    e <- GiottoClass::createNetwork(locs, GiottoClass::kNNNetworkParam(
        k = k, include_weight = FALSE, output = "data.table"))
    e[, weight := exp(-distance^2 / stats::median(distance)^2), by = from]
    e[!is.finite(weight), weight := 1] # all neighbours at distance 0
    e[, weight := weight / sum(weight), by = from]
    i <- match(e$from, ids)
    j <- match(e$to, ids)
    phi <- atan2(locs[j, 2] - locs[i, 2], locs[j, 1] - locs[i, 1])
    n <- length(ids)
    sp <- function(v) Matrix::sparseMatrix(i = i, j = j, x = v,
        dims = c(n, n), dimnames = list(ids, ids))
    list(w = sp(e$weight), wc = sp(e$weight * cos(phi)),
        ws = sp(e$weight * sin(phi)), a = sp(rep(1 / k, length(i))))
}

.neighbor_weight_matrix <- function(locs, k) .neighbor_operators(locs, k)$w

# |sum_i w_ji exp(i phi_ji) (x_i - mean of j's neighbours)| per feature, as
# Banksy's harmonic m = 1 with centring; features x cells
.neighbor_agf <- function(x, ops) {
    m <- as.matrix(x %*% Matrix::t(ops$a))
    re <- as.matrix(x %*% Matrix::t(ops$wc)) -
        sweep(m, 2L, Matrix::rowSums(ops$wc), "*")
    im <- as.matrix(x %*% Matrix::t(ops$ws)) -
        sweep(m, 2L, Matrix::rowSums(ops$ws), "*")
    sqrt(re^2 + im^2)
}

# z-score each row (feature) across columns (cells), sample SD; constant rows
# become 0 (Banksy's scaler).
.scale_rows <- function(x) {
    mu <- rowMeans(x)
    s <- MatrixGenerics::rowSds(x)
    out <- (x - mu) / s
    out[!is.finite(out)] <- 0
    out
}


# runNeighborhoodPCA ####

#' @title Neighborhood PCA (BANKSY embedding)
#' @name runNeighborhoodPCA
#' @description
#' Computes the BANKSY embedding of Singhal et al. (2024) natively and stores
#' it as a dimension reduction: each cell's values are joined by the weighted
#' mean of its spatial neighbours, weighted by `lambda`, and reduced by PCA
#' (see [neighborhoodPCAParam()]). Cluster the result with
#' [createNearestNetwork()] and [doLeidenCluster()] (`lambda` about 0.8 for
#' spatial domains, about 0.2 for neighbourhood-aware cell types).
#'
#' This reimplements BANKSY's embedding, including its gradient (AGF)
#' features with `compute_agf = TRUE`. **If you use it, please cite the
#' original paper** (see References). `GiottoConnect::runBanksy()` runs the original \pkg{Banksy}
#' package.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param input `"expression"` (default) or `"dim_reduction"`: the values the
#'   neighbourhood is computed on
#' @param expression_values expression values to use with
#'   `input = "expression"`
#' @param feats_to_use features to use with `input = "expression"`: a
#'   character vector, or the name of a logical column of the feature
#'   metadata (e.g. `"hvf"`). `NULL` uses all features.
#' @param dim_reduction_to_use,dim_reduction_name,dimensions_to_use the
#'   embedding to use with `input = "dim_reduction"`
#' @param spat_loc_name name of the spatial locations
#' @param lambda,k_geom,npcs,compute_agf,absolute,seed see
#'   [neighborhoodPCAParam()]
#' @param name name of the stored dimension reduction (reduction method
#'   `"neighborhood_pca"`)
#' @param verbose be verbose
#' @returns The giotto object with the embedding stored as a dimension
#'   reduction (method `"neighborhood_pca"`, name `name`).
#' @references Singhal V, Chou N, Lee J, et al. BANKSY unifies cell typing
#'   and tissue domain segmentation for scalable spatial omics data analysis.
#'   *Nature Genetics* 56, 431–441 (2024). \doi{10.1038/s41588-024-01664-3}
#' @seealso [neighborhoodPCAParam()], `GiottoConnect::runBanksy()`
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- runNeighborhoodPCA(g, feats_to_use = "hvf", lambda = 0.8)
#' g <- createNearestNetwork(g, dim_reduction_to_use = "neighborhood_pca",
#'     dim_reduction_name = "neighborhood_pca", name = "neighborhood_pca")
#' g <- doLeidenCluster(g, network_name = "neighborhood_pca",
#'     name = "nbr_domain")
#' @md
#' @export
runNeighborhoodPCA <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        input = c("expression", "dim_reduction"),
        expression_values = "normalized",
        feats_to_use = NULL,
        dim_reduction_to_use = "pca",
        dim_reduction_name = NULL,
        dimensions_to_use = NULL,
        spat_loc_name = "raw",
        lambda = 0.2,
        k_geom = 15L,
        npcs = 20L,
        compute_agf = FALSE,
        absolute = FALSE,
        seed = 1234L,
        name = "neighborhood_pca",
        verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    input <- match.arg(input)
    checkmate::assert_string(name)
    param <- neighborhoodPCAParam(lambda = lambda, k_geom = k_geom,
        npcs = npcs, compute_agf = compute_agf, absolute = absolute,
        seed = seed)
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    x <- switch(input,
        expression = {
            ex <- getExpression(gobject, spat_unit = spat_unit,
                feat_type = feat_type, values = expression_values,
                output = "exprObj")
            feats <- .npca_feats(gobject, feats_to_use, spat_unit, feat_type)
            .read_expr_matrix(ex[], feats)
        },
        dim_reduction = {
            emb <- getDimReduction(gobject, spat_unit = spat_unit,
                feat_type = feat_type, reduction = "cells",
                reduction_method = dim_reduction_to_use,
                name = dim_reduction_name, output = "matrix")
            if (!is.null(dimensions_to_use)) {
                emb <- emb[, intersect(seq_len(ncol(emb)), dimensions_to_use),
                    drop = FALSE]
            }
            t(emb)
        }
    )
    sl <- getSpatialLocations(gobject, spat_unit = spat_unit,
        name = spat_loc_name, output = "data.table")
    locs <- as.matrix(sl[, c("sdimx", "sdimy")])
    rownames(locs) <- sl$cell_ID

    vmsg(.v = verbose, "runNeighborhoodPCA: neighbourhood and PCA")
    res <- reduceData(x, param, spatial_locs = locs)
    coords <- res$u
    dimnames(coords) <- list(colnames(x), paste0("Dim.", seq_len(ncol(coords))))

    dim_obj <- createDimObj(
        coordinates = coords, name = name, spat_unit = spat_unit,
        feat_type = feat_type, method = "neighborhood_pca",
        reduction = "cells",
        misc = list(input = input, eigenvalues = res$eigenvalues,
            loadings = res$v)
    )
    gobject <- setDimReduction(gobject, x = dim_obj, spat_unit = spat_unit,
        feat_type = feat_type, reduction = "cells",
        reduction_method = "neighborhood_pca", name = name, verbose = FALSE)
    update_giotto_params(gobject, description = "_neighborhood_pca",
        attachments = list(param = param))
}

# features to use: NULL (all), a character vector, or a logical feature
# metadata column
.npca_feats <- function(gobject, feats_to_use, spat_unit, feat_type) {
    if (is.null(feats_to_use)) return(NULL)
    if (length(feats_to_use) == 1L) {
        fm <- fDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
        if (feats_to_use %in% colnames(fm)) {
            col <- fm[[feats_to_use]]
            keep <- if (is.logical(col)) col else col %in% c("yes", "TRUE")
            return(fm$feat_ID[keep])
        }
    }
    feats_to_use
}
