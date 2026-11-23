#' @include clustering.R
NULL

# hmrfParam ####
#
# HMRF spatial domains as a `clusterData()` verb (adr/0008), the sibling of
# `deconvParam` in deconv-param.R. The algorithm is not in this package: it is
# the spatialHMRF package (github: drieslab/spatialHMRF), which Giotto
# suggests, the way deconvolution uses spatialDWLS. Giotto owns the param
# class, the verb and the gobject plumbing in runHMRF() and
# findSpatialDomains(); the `clusterData` method below is the only code here
# that calls the engine.

#' @title HMRF clustering parameters
#' @name hmrfParam
#' @description
#' Parameters for [clusterData()] that fit the hidden Markov random field
#' (HMRF) model of Zhu et al. (2018), as in the \pkg{smfishHmrf} package used
#' by [doHMRF_V2()]: each cell's z-scored expression of spatial genes follows
#' a multivariate Gaussian per domain, and a Potts prior of strength `beta`
#' favours domains shared with spatial neighbours. Fitted by EM with
#' iterated conditional modes over the colour classes of the spatial network,
#' from a k-means start.
#'
#' The model is fitted by the \pkg{spatialHMRF} package, a reimplementation
#' of \pkg{smfishHmrf}. **If you use it, please cite the original paper**
#' (see References).
#' @param k integer, number of domains
#' @param beta numeric, strength of the spatial prior (0: no prior)
#' @param zscore `"cell"` (default, as [initHMRF_V2()]): z-score each cell
#'   across the features after the input values; `"none"`: use the values as
#'   given
#' @param nstart integer, k-means starts
#' @param seed integer, seed of the k-means start
#' @param maxit integer, maximum EM iterations
#' @param err numeric, relative change of means and covariances at which EM
#'   stops
#' @param n_icm integer, ICM sweeps per EM iteration
#' @param damp_tolerance,damp_factor a domain's covariance whose determinant
#'   is below `damp_tolerance` gets the smallest diagonal term (searched in
#'   steps of `damp_factor` from 1e-4) that brings it above, as in
#'   \pkg{smfishHmrf}
#' @returns An `hmrfParam` object.
#' @references Zhu Q, Shah S, Dries R, Cai L, Yuan GC. Identification of
#'   spatially associated subpopulations by combining scRNAseq and sequential
#'   fluorescence in situ hybridization data. *Nature Biotechnology* 36,
#'   1183–1190 (2018). \doi{10.1038/nbt.4260}
#' @seealso [runHMRF()], [doHMRF_V2()]
#' @examples
#' p <- hmrfParam(k = 5, beta = 2)
#' @md
#' @export
hmrfParam <- function(k = 10L, beta = 2, zscore = c("cell", "none"),
    nstart = 1000L, seed = 100L, maxit = 50L, err = 1e-7, n_icm = 10L,
    damp_tolerance = 1e-5, damp_factor = 1.05) {
    checkmate::assert_int(k, lower = 2)
    checkmate::assert_number(beta, lower = 0)
    zscore <- match.arg(zscore)
    checkmate::assert_count(nstart, positive = TRUE)
    checkmate::assert_count(seed)
    checkmate::assert_count(maxit, positive = TRUE)
    checkmate::assert_number(err, lower = 0)
    checkmate::assert_count(n_icm, positive = TRUE)
    new("hmrfParam", param = list(k = as.integer(k), beta = beta,
        zscore = zscore, nstart = as.integer(nstart), seed = as.integer(seed),
        maxit = as.integer(maxit), err = err, n_icm = as.integer(n_icm),
        damp_tolerance = damp_tolerance, damp_factor = damp_factor))
}

#' @rdname hmrfParam
#' @exportClass hmrfParam
setClass("hmrfParam", contains = "BlusterParam", slots = list(param = "list"))

#' @rdname hmrfParam
#' @param x cells x features values (`matrix` or `Matrix`), cell IDs as row
#'   names (the orientation of the other `clusterData()` methods)
#' @param param an `hmrfParam`
#' @param ... unused
#' @param network spatial network: a `data.frame` with columns `from` and `to`
#'   (cell IDs), each pair listed once
#' @param blocks optional list of integer vectors (positions in
#'   `rownames(x)`) that partition the cells into independent sets of the
#'   network; a DSatur colouring of the network when `NULL`
#' @param init optional list with `mu` (features x k), `sigma`
#'   (features x features x k) and `damp` (length k), e.g. to reuse one start
#'   for several `beta`
#' @details The computation is done by the \pkg{spatialHMRF} package
#'   (`spatialHMRF::hmrf()`), which is suggested rather than imported:
#'   install it with `remotes::install_github("drieslab/spatialHMRF")`.
#'   `clusterData(x, param, network = )` returns a `data.table` with
#'   `cell_ID` and `cluster` (integer). Its attributes: `prob` (cells x k
#'   posterior probabilities), `mu`, `sigma`, `likelihood`, `iterations`,
#'   `init` and `blocks`.
#' @export
setMethod("clusterData", signature("allMatrix", "hmrfParam"),
    function(x, param, ..., network, blocks = NULL, init = NULL) {
        .hmrf_engine_check()
        p <- param@param
        fit <- spatialHMRF::hmrf(as.matrix(x), network, k = p$k,
            beta = p$beta, zscore = p$zscore, init = init, blocks = blocks,
            nstart = p$nstart, seed = p$seed, maxit = p$maxit, err = p$err,
            n_icm = p$n_icm, damp_tolerance = p$damp_tolerance,
            damp_factor = p$damp_factor)
        res <- data.table::data.table(cell_ID = rownames(x),
            cluster = fit$class)
        for (a in c("prob", "mu", "sigma", "likelihood", "iterations",
            "init", "blocks")) {
            attr(res, a) <- fit[[a]]
        }
        res
    }
)

# The engine is suggested, not imported: stop with the install line when it
# is missing.
#' @keywords internal
#' @noRd
.hmrf_engine_check <- function() {
    package_check(
        pkg_name = "spatialHMRF",
        repository = "github:drieslab/spatialHMRF"
    )
}


# runHMRF ####

#' @title Spatial domains with HMRF
#' @name runHMRF
#' @description
#' Finds spatial domains with the hidden Markov random field model of Zhu et
#' al. (2018) (see [hmrfParam()]), natively, for one or more `beta`. For each
#' `beta` the domains are added to the cell metadata (column
#' `<name>_k<k>_b<beta>`) and the posterior domain probabilities stored as a
#' spatial enrichment of the same name. All `beta` share one k-means start.
#'
#' This reimplements the model of \pkg{smfishHmrf} ([doHMRF_V2()]) without its
#' spatial gene selection: pass the features with `feats_to_use`. **If you use
#' it, please cite the original paper** (see References).
#' @param gobject giotto object
#' @param spat_unit,feat_type spatial unit and feature type
#' @param expression_values expression values (default `"scaled"`, as
#'   [initHMRF_V2()])
#' @param feats_to_use features: a character vector (e.g. spatial genes from
#'   [binSpect()]) or the name of a logical feature metadata column such as
#'   `"hvf"`
#' @param spatial_network_name spatial network
#' @param k,zscore,nstart,seed,maxit,err,n_icm see [hmrfParam()]
#' @param betas numeric vector of `beta` values
#' @param name prefix of the results
#' @param verbose be verbose
#' @returns The giotto object with the domains and probabilities.
#' @references Zhu Q, Shah S, Dries R, Cai L, Yuan GC. Identification of
#'   spatially associated subpopulations by combining scRNAseq and sequential
#'   fluorescence in situ hybridization data. *Nature Biotechnology* 36,
#'   1183–1190 (2018). \doi{10.1038/nbt.4260}
#' @seealso [hmrfParam()], [doHMRF_V2()]
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- runHMRF(g, feats_to_use = "hvf", k = 5, betas = c(0, 2),
#'     nstart = 20)
#' @md
#' @export
runHMRF <- function(gobject, spat_unit = NULL, feat_type = NULL,
    expression_values = "scaled", feats_to_use = NULL,
    spatial_network_name = "Delaunay_network", k = 10L, betas = 2,
    zscore = c("cell", "none"), nstart = 1000L, seed = 100L, maxit = 50L,
    err = 1e-7, n_icm = 10L, name = "hmrf", verbose = NULL) {
    cell_ID <- NULL # NSE
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_numeric(betas, lower = 0, min.len = 1L)
    zscore <- match.arg(zscore)
    spat_unit <- set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject = gobject,
        spat_unit = spat_unit, feat_type = feat_type)

    ex <- getExpression(gobject, spat_unit = spat_unit, feat_type = feat_type,
        values = expression_values, output = "exprObj")
    feats <- .npca_feats(gobject, feats_to_use, spat_unit, feat_type)
    x <- t(as.matrix(.read_expr_matrix(ex[], feats)))
    net <- getSpatialNetwork(gobject, spat_unit = spat_unit,
        name = spatial_network_name, output = "networkDT")

    init <- NULL
    blocks <- NULL
    for (b in betas) {
        param <- hmrfParam(k = k, beta = b, zscore = zscore, nstart = nstart,
            seed = seed, maxit = maxit, err = err, n_icm = n_icm)
        vmsg(.v = verbose, sprintf("runHMRF: k = %d, beta = %g", k, b))
        res <- clusterData(x, param, network = net, blocks = blocks,
            init = init)
        init <- attr(res, "init")
        blocks <- attr(res, "blocks")
        key <- sprintf("%s_k%d_b%s", name, k, format(b))

        new_meta <- data.table::data.table(cell_ID = res$cell_ID)
        new_meta[[key]] <- res$cluster
        gobject <- addCellMetadata(gobject, spat_unit = spat_unit,
            feat_type = feat_type, new_metadata = new_meta,
            by_column = TRUE, column_cell_ID = "cell_ID")
        prob <- data.table::as.data.table(attr(res, "prob"))
        data.table::setnames(prob, paste0("domain_", seq_len(k)))
        prob[, cell_ID := res$cell_ID]
        enr <- createSpatEnrObj(prob, name = key, spat_unit = spat_unit,
            feat_type = feat_type, method = "HMRF",
            misc = list(likelihood = attr(res, "likelihood"),
                iterations = attr(res, "iterations"), beta = b),
            verbose = FALSE)
        gobject <- setSpatialEnrichment(gobject, x = enr,
            spat_unit = spat_unit, feat_type = feat_type, name = key,
            verbose = FALSE)
        gobject <- update_giotto_params(gobject,
            description = paste0("_hmrf_b", format(b)),
            attachments = list(param = param))
    }
    gobject
}


# findSpatialDomains ####

#' @title Find spatial domains
#' @name findSpatialDomains
#' @description
#' One command for spatial domains, with two native methods:
#' * `"neighborhood"`: the BANKSY embedding ([runNeighborhoodPCA()], Singhal
#'   et al. 2024) at `lambda`, clustered by Leiden on an sNN network;
#' * `"hmrf"`: the hidden Markov random field ([runHMRF()], Zhu et al. 2018).
#'
#' **If you use them, please cite the original papers** (see the pages of
#' [runNeighborhoodPCA()] and [runHMRF()]).
#' @param gobject giotto object
#' @param method `"neighborhood"` or `"hmrf"`
#' @param spat_unit,feat_type spatial unit and feature type
#' @param feats_to_use features (character vector or a logical feature
#'   metadata column such as `"hvf"`)
#' @param k number of domains (`"hmrf"`)
#' @param beta spatial prior strength (`"hmrf"`)
#' @param lambda neighbourhood weight (`"neighborhood"`)
#' @param resolution Leiden resolution (`"neighborhood"`)
#' @param name name of the cell metadata column with the domains
#' @param ... passed to [runNeighborhoodPCA()] or [runHMRF()] (named
#'   arguments of [neighborhoodPCAParam()] or [hmrfParam()])
#' @returns The giotto object with the domains in cell metadata `name`.
#' @references Singhal V, Chou N, Lee J, et al. BANKSY unifies cell typing
#'   and tissue domain segmentation for scalable spatial omics data analysis.
#'   *Nature Genetics* 56, 431–441 (2024). \doi{10.1038/s41588-024-01664-3}
#'
#'   Zhu Q, Shah S, Dries R, Cai L, Yuan GC. Identification of spatially
#'   associated subpopulations by combining scRNAseq and sequential
#'   fluorescence in situ hybridization data. *Nature Biotechnology* 36,
#'   1183–1190 (2018). \doi{10.1038/nbt.4260}
#' @md
#' @export
findSpatialDomains <- function(gobject, method = c("neighborhood", "hmrf"),
    spat_unit = NULL, feat_type = NULL, feats_to_use = "hvf", k = 10L,
    beta = 2, lambda = 0.8, resolution = 0.8, name = "spatial_domain", ...) {
    method <- match.arg(method)
    spat_unit <- set_default_spat_unit(gobject = gobject, spat_unit = spat_unit)
    feat_type <- set_default_feat_type(gobject = gobject,
        spat_unit = spat_unit, feat_type = feat_type)
    # inner wrappers do not record history (nested calls); one entry below
    quiet <- function(code) gwith_options(list(giotto.update_param = FALSE),
        code)
    if (method == "neighborhood") {
        param <- neighborhoodPCAParam(lambda = lambda, ...)
        gobject <- quiet(runNeighborhoodPCA(gobject, spat_unit = spat_unit,
            feat_type = feat_type, feats_to_use = feats_to_use,
            lambda = lambda, name = name, ...))
        gobject <- createNearestNetwork(gobject, spat_unit = spat_unit,
            feat_type = feat_type, type = "sNN",
            dim_reduction_to_use = "neighborhood_pca",
            dim_reduction_name = name, name = name, verbose = FALSE)
        gobject <- quiet(doLeidenCluster(gobject, spat_unit = spat_unit,
            feat_type = feat_type, nn_network_to_use = "sNN",
            network_name = name, name = name, resolution = resolution,
            n_iterations = -1))
        return(update_giotto_params(gobject, description = "_spatial_domains",
            attachments = list(param = param)))
    }
    param <- hmrfParam(k = k, beta = beta, ...)
    gobject <- quiet(runHMRF(gobject, spat_unit = spat_unit,
        feat_type = feat_type, feats_to_use = feats_to_use, k = k,
        betas = beta, name = name, ...))
    key <- sprintf("%s_k%d_b%s", name, k, format(beta))
    meta <- pDataDT(gobject, spat_unit = spat_unit, feat_type = feat_type)
    new_meta <- data.table::data.table(cell_ID = meta$cell_ID,
        domain = meta[[key]])
    data.table::setnames(new_meta, "domain", name)
    gobject <- addCellMetadata(gobject, spat_unit = spat_unit,
        feat_type = feat_type, new_metadata = new_meta, by_column = TRUE,
        column_cell_ID = "cell_ID")
    update_giotto_params(gobject, description = "_spatial_domains",
        attachments = list(param = param))
}
