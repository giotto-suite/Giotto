#' @include clustering.R
NULL

# hmrfParam ####
#
# Native hidden Markov random field (HMRF) clustering, the model of Zhu et al.
# 2018 as implemented in smfishHmrf (doHMRF_V2), as a `clusterData()` verb
# (adr/0008): multivariate Gaussian emissions on z-scored spatial genes, a
# Potts prior over the spatial network with strength beta, k-means
# initialisation and EM with iterated conditional modes (ICM) over the colour
# classes of the network. Each colour class is an independent set, so its ICM
# update and the neighbour counts are sparse matrix products, which give
# smfishHmrf's result without its per-cell loops.

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
#' This is a reimplementation of that model. **If you use it, please cite
#' the original paper** (see References).
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
#' @details `clusterData(x, param, network = )` returns a `data.table` with
#'   `cell_ID` and `cluster` (integer). Its attributes: `prob` (cells x k
#'   posterior probabilities), `mu`, `sigma`, `likelihood`, `iterations` and
#'   `init`.
#' @export
setMethod("clusterData", signature("allMatrix", "hmrfParam"),
    function(x, param, ..., network, blocks = NULL, init = NULL) {
        p <- param@param
        ids <- rownames(x)
        y <- as.matrix(x)
        if (identical(p$zscore, "cell")) y <- t(scale(t(y)))
        y[!is.finite(y)] <- 0
        adj <- .hmrf_adjacency(network, ids)
        if (is.null(blocks)) blocks <- .hmrf_blocks(adj)
        if (is.null(init)) init <- .hmrf_init(y, p)
        fit <- .hmrf_fit(y, adj, blocks, init, p)
        res <- data.table::data.table(cell_ID = ids, cluster = fit$class)
        attr(res, "prob") <- `rownames<-`(fit$prob, ids)
        for (a in c("mu", "sigma", "likelihood", "iterations")) {
            attr(res, a) <- fit[[a]]
        }
        attr(res, "init") <- init
        res
    }
)

# symmetric neighbour counts from an edge list of cell IDs
.hmrf_adjacency <- function(network, ids) {
    i <- match(as.character(network$from), ids)
    j <- match(as.character(network$to), ids)
    keep <- !is.na(i) & !is.na(j) & i != j
    i <- i[keep]
    j <- j[keep]
    n <- length(ids)
    Matrix::sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(n, n))
}

# independent sets: colour classes of a DSatur colouring
.hmrf_blocks <- function(adj) {
    g <- igraph::graph_from_adjacency_matrix(adj > 0, mode = "undirected",
        diag = FALSE)
    col <- igraph::greedy_vertex_coloring(g, heuristic = "dsatur")
    lapply(sort(unique(col)), function(cl) which(col == cl))
}

# k-means start, per-domain covariances and damping (smfishHmrf defaults)
.hmrf_init <- function(y, p) {
    if (p$seed > 0) set.seed(p$seed)
    kk <- stats::kmeans(y, p$k, nstart = p$nstart, iter.max = 100)
    m <- ncol(y)
    sigma <- array(0, c(m, m, p$k))
    damp <- numeric(p$k)
    for (j in seq_len(p$k)) {
        sigma[, , j] <- stats::cov(y[kk$cluster == j, , drop = FALSE])
        d <- .hmrf_damp(sigma[, , j], p$damp_factor, p$damp_tolerance)
        damp[j] <- if (is.null(d)) 0 else d
    }
    list(mu = t(kk$centers), sigma = sigma, damp = damp)
}

.hmrf_damp <- function(sigma, factor, cutoff, start = 1e-4) {
    if (det(sigma) >= cutoff) return(NULL)
    eye <- diag(nrow(sigma))
    damp <- start
    if (det(sigma + eye * damp) < cutoff) {
        repeat {
            damp <- damp * factor
            if (det(sigma + eye * damp) > cutoff) break
        }
    } else {
        repeat {
            damp <- damp / factor
            if (det(sigma + eye * damp) < cutoff) {
                damp <- damp * factor
                break
            }
        }
    }
    damp
}

# Gaussian densities of every cell under every domain (clipped as smfishHmrf)
.hmrf_density <- function(y, mu, sigma, damp) {
    n <- nrow(y)
    m <- ncol(y)
    k <- ncol(mu)
    den <- matrix(0, n, k)
    for (j in seq_len(k)) {
        s <- sigma[, , j] + diag(m) * damp[j]
        dist <- sweep(y, 2L, mu[, j])
        e <- -0.5 * rowSums((dist %*% solve(s)) * dist) -
            ((m / 2) * log(2 * pi) + 0.5 * log(det(s)))
        den[, j] <- exp(pmin(700, pmax(e, -700)))
    }
    den
}

.hmrf_onehot <- function(cl, k) {
    z <- matrix(0, length(cl), k)
    z[cbind(seq_along(cl), cl)] <- 1
    z
}

.hmrf_relerr <- function(a, b) max(abs(a - b)) / (1 + max(abs(a), abs(b)))

.hmrf_fit <- function(y, adj, blocks, init, p) {
    k <- p$k
    n <- nrow(y)
    mu <- init$mu
    sigma <- init$sigma
    damp <- init$damp
    z <- .hmrf_onehot(max.col(.hmrf_density(y, mu, sigma, damp)), k)
    it <- 0L
    repeat {
        it <- it + 1L
        mu_old <- mu
        sigma_old <- sigma
        den <- .hmrf_density(y, mu, sigma, damp)
        logden <- log(den)
        for (s in seq_len(p$n_icm)) {
            for (b in blocks) {
                ni <- as.matrix(adj[b, , drop = FALSE] %*% z)
                cl <- max.col(logden[b, , drop = FALSE] + p$beta * ni,
                    ties.method = "first")
                z[b, ] <- .hmrf_onehot(cl, k)
            }
        }
        nij <- as.matrix(adj %*% z)
        prob <- den * exp(p$beta * nij)
        prob <- prob / rowSums(prob)
        bad <- !is.finite(rowSums(prob))
        prob[bad, ] <- 1 / k
        unnorm <- den * exp(pmin(700, pmax(p$beta * nij, -700)))
        w <- colSums(prob)
        mu <- crossprod(y, prob) %*% diag(1 / w, k)
        for (j in seq_len(k)) {
            dist <- sweep(y, 2L, mu[, j])
            sigma[, , j] <- crossprod(dist * prob[, j], dist) / w[j]
        }
        if ((.hmrf_relerr(mu_old, mu) < p$err &&
            .hmrf_relerr(sigma_old, sigma) < p$err) || it > p$maxit) break
    }
    t_prob <- colMeans(unnorm / rowSums(unnorm))
    list(class = max.col(prob, "first"), prob = prob, mu = mu, sigma = sigma,
        likelihood = sum(log(unnorm %*% t_prob)), iterations = it)
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
        if (is.null(blocks)) {
            blocks <- .hmrf_blocks(.hmrf_adjacency(net, rownames(x)))
        }
        res <- clusterData(x, param, network = net, blocks = blocks,
            init = init)
        init <- attr(res, "init")
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
