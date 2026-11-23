#' Multi omics integration with block-weighted PCs
#'
#' @description
#' Joint embedding of two or more feature types from their own dimension
#' reductions (multiblock PCA with block scaling; Westerhuis et al. 1998,
#' *J Chemometrics* 12:301). Each feature type's coordinates (a block) are
#' divided by the square root of their total variance, so no block dominates
#' by its number of dimensions or its scale, multiplied by the square root of
#' its weight and concatenated. With `rotate = TRUE` the concatenation is
#' rotated to its principal axes, so its dimensions are ordered by variance as
#' in a PCA; a rotation keeps all distances, so neighbors computed on all
#' dimensions are the same either way.
#'
#' Memory and time are linear in the number of cells (one N x d matrix, d the
#' dimensions summed over the feature types), so this scales to whole slides
#' and works on GiottoDisk-backed objects.
#'
#' @param gobject A Giotto object with a dimension reduction (e.g. PCA) for
#'   each feature type in `feat_types`.
#' @param spat_unit spatial unit
#' @param feat_types feature types to integrate (two or more). The joint
#'   embedding is stored under each of them.
#' @param reduction_methods reduction method of each feature type's embedding.
#'   Default `c("pca", "pca")`.
#' @param reduction_names name of each feature type's embedding. `NULL`
#'   (default) uses the names `runPCA()` gives: `"pca"` for `rna`,
#'   `"<feat_type>.pca"` otherwise. `list_dim_reductions(gobject)` lists the
#'   stored names.
#' @param dimensions_to_use dimensions (indices) of each embedding to use: one
#'   integer vector for all feature types, or a list with one per feature type
#'   (named by feature type, or in the order of `feat_types`). Indices beyond
#'   those stored are dropped. `NULL` (default): all stored dimensions.
#' @param block_weights share of the joint variance given to each feature
#'   type. `NULL` (default) weights them equally; otherwise one positive number
#'   per feature type (named by feature type, or in the order of
#'   `feat_types`), rescaled to sum to 1.
#' @param rotate logical. Rotate the joint embedding to its principal axes
#'   (default `TRUE`).
#' @param name name of the stored embedding. Default `"blockpca"`.
#' @param integrated_k number of neighbors of each cell in the integrated kNN,
#'   searched on all dimensions of the joint embedding. Default 20.
#' @param integrated_network_name name of the integrated kNN network, stored as
#'   a `kNN` network under `feat_types[1]`. Default `"<name>_kNN"`.
#' @param integrated_snn shared nearest neighbor network of the integrated kNN
#'   to store for clustering, as an `sNN` network named `integrated_snn_name`
#'   under `feat_types[1]`, as in [runWNN()]. `"none"` (default) stores no sNN;
#'   `"giotto"` keeps the kNN edges whose cells share at least 5 neighbors,
#'   plus each cell's 3 strongest (as `createNearestNetwork(type = "sNN")`);
#'   `"jaccard"` weights every kNN edge by the Jaccard overlap of the two
#'   cells' neighbor sets and drops edges below `snn_prune`.
#' @param snn_prune smallest Jaccard overlap kept with
#'   `integrated_snn = "jaccard"`. Default 1/15.
#' @param integrated_snn_name name of the sNN network. `NULL` (default):
#'   `integrated_network_name` with `"kNN"` replaced by `"sNN"`.
#' @param engine neighbor search: `"dbscan"` (exact, default) or `"hnsw"`
#'   (approximate, faster on large data).
#' @param n_threads threads of the `"hnsw"` search. `NULL`: its default.
#' @param verbose be verbose
#'
#' @details
#' For block b with selected coordinates X_b (cells x d_b) and weight w_b,
#' the joint embedding is `[X_1 s_1, ..., X_B s_B]` with
#' `s_b = sqrt(w_b / sum_j var(X_b[, j]))`, so block b holds w_b of the total
#' variance (1). The rotation centers the joint embedding and projects it on
#' the eigenvectors of its covariance, each signed so that its largest entry
#' is positive.
#'
#' Stored in the `misc` slot of each dimension reduction: `eigenvalues` (the
#' variance of each joint dimension; `NULL` without rotation) and
#' `block_input`, with the feature types, embeddings and dimensions used, the
#' weights, the scale factors s_b, and `block_share`, the share of each joint
#' dimension's variance that comes from each block.
#'
#' Unlike [runWNN()], the weights are the same for every cell.
#'
#' @returns The Giotto object with the joint embedding stored, for every
#'   feature type in `feat_types`, as reduction method `"blockpca"` named
#'   `name`; the integrated kNN (`integrated_network_name`) and, unless
#'   `integrated_snn = "none"`, its sNN (`integrated_snn_name`) under
#'   `feat_types[1]`. Cluster with `doLeidenCluster(feat_type = feat_types[1],
#'   nn_network_to_use = "sNN", network_name = integrated_snn_name)` and embed
#'   with `runIntegratedUMAP(graph = "knn", network_name =
#'   integrated_network_name)`, as for [runWNN()].
#' @seealso [runWNN()], [runIntegratedUMAP()]
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' x <- getExpression(g, values = "raw", output = "matrix")[1:30, ]
#' g <- setExpression(g, createExprObj(x, name = "raw", feat_type = "protein"))
#' g <- normalizeGiotto(g, feat_type = "protein", verbose = FALSE)
#' g <- runPCA(g, feat_type = "protein", feats_to_use = NULL, ncp = 10)
#'
#' g <- runBlockPCA(g, feat_types = c("rna", "protein"),
#'     dimensions_to_use = list(1:10, 1:10), integrated_k = 10,
#'     integrated_snn = "giotto")
#' getDimReduction(g, reduction_method = "blockpca", name = "blockpca",
#'     output = "matrix")[1:3, 1:3]
#' g <- doLeidenCluster(g, nn_network_to_use = "sNN",
#'     network_name = "blockpca_sNN", name = "cluster_block")
#' g <- runIntegratedUMAP(g, graph = "knn", network_name = "blockpca_kNN",
#'     k = 10, name = "umap_block")
#' @md
#' @export
runBlockPCA <- function(
        gobject,
        spat_unit = "cell",
        feat_types = c("rna", "protein"),
        reduction_methods = c("pca", "pca"),
        reduction_names = NULL,
        dimensions_to_use = NULL,
        block_weights = NULL,
        rotate = TRUE,
        name = "blockpca",
        integrated_k = 20,
        integrated_network_name = paste0(name, "_kNN"),
        integrated_snn = c("none", "giotto", "jaccard"),
        snn_prune = 1 / 15,
        integrated_snn_name = NULL,
        engine = c("dbscan", "hnsw"),
        n_threads = NULL,
        verbose = FALSE) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_string(spat_unit)
    checkmate::assert_character(feat_types, min.len = 2L, unique = TRUE,
        any.missing = FALSE)
    n_mod <- length(feat_types)
    reduction_methods <- .wnn_per_mod(reduction_methods, feat_types,
        "reduction_methods", fun = "runBlockPCA")
    if (is.null(reduction_names)) {
        reduction_names <- ifelse(feat_types == "rna", reduction_methods,
            paste0(feat_types, ".", reduction_methods))
    }
    reduction_names <- .wnn_per_mod(reduction_names, feat_types,
        "reduction_names", fun = "runBlockPCA")
    dimensions_to_use <- .wnn_per_mod_dims(dimensions_to_use, feat_types,
        fun = "runBlockPCA")
    block_weights <- .block_weights(block_weights, feat_types)
    checkmate::assert_flag(rotate)
    checkmate::assert_string(name)
    checkmate::assert_count(integrated_k, positive = TRUE)
    checkmate::assert_string(integrated_network_name)
    integrated_snn <- match.arg(integrated_snn)
    checkmate::assert_number(snn_prune, lower = 0, upper = 1)
    if (is.null(integrated_snn_name)) {
        integrated_snn_name <- sub("kNN", "sNN", integrated_network_name,
            fixed = TRUE)
    }
    checkmate::assert_string(integrated_snn_name)
    engine <- match.arg(engine)

    .wnn_check_reductions(gobject, spat_unit, feat_types, "runBlockPCA")

    blocks <- dims <- list()
    scale <- numeric(n_mod)
    cells <- NULL
    for (i in seq_len(n_mod)) {
        ft <- feat_types[i]
        emb <- getDimReduction(gobject,
            spat_unit = spat_unit, feat_type = ft, reduction = "cells",
            reduction_method = reduction_methods[i],
            name = reduction_names[i], output = "matrix"
        )
        if (is.null(rownames(emb)) || anyDuplicated(rownames(emb))) {
            stop("[runBlockPCA] the embedding of ", ft, " needs unique cell ",
                "IDs as row names", call. = FALSE)
        }
        all_d <- seq_len(ncol(emb))
        d <- dimensions_to_use[[i]]
        d <- if (is.null(d)) all_d else intersect(as.integer(d), all_d)
        if (!length(d)) {
            stop("[runBlockPCA] no stored dimensions of ", ft, " are selected",
                call. = FALSE)
        }
        vmsg(.v = verbose, sprintf("runBlockPCA: %s uses dimensions %s", ft,
            .wnn_range_txt(d)))
        emb <- emb[, d, drop = FALSE]
        storage.mode(emb) <- "double"
        bad <- rowSums(!is.finite(emb)) > 0L
        if (any(bad)) {
            stop(sprintf(paste0("[runBlockPCA] the embedding of %s has NA ",
                "or infinite values for %d cells"), ft, sum(bad)),
                call. = FALSE)
        }

        if (is.null(cells)) {
            cells <- rownames(emb)
        } else if (!setequal(cells, rownames(emb))) {
            n_miss <- sum(!cells %in% rownames(emb))
            n_extra <- sum(!rownames(emb) %in% cells)
            stop(sprintf(paste0("[runBlockPCA] %s and %s embed different ",
                "cells (%d missing from %s, %d only in %s); subset them to ",
                "the same cells first"), feat_types[1], ft, n_miss, ft,
                n_extra, ft), call. = FALSE)
        } else if (!identical(rownames(emb), cells)) {
            emb <- emb[cells, , drop = FALSE]
        }

        total_var <- sum(apply(emb, 2L, stats::var))
        if (!is.finite(total_var) || total_var <= 0) {
            stop("[runBlockPCA] the embedding of ", ft, " has no variance",
                call. = FALSE)
        }
        scale[i] <- sqrt(block_weights[i] / total_var)
        blocks[[ft]] <- emb * scale[i]
        dims[[ft]] <- d
    }

    joint <- do.call(cbind, unname(blocks))
    block_of <- rep(feat_types, vapply(blocks, ncol, integer(1L)))
    eigenvalues <- NULL
    loading <- diag(ncol(joint))
    if (rotate) {
        vmsg(.v = verbose, "runBlockPCA: rotating to principal axes")
        joint <- sweep(joint, 2L, colMeans(joint))
        e <- eigen(crossprod(joint) / (nrow(joint) - 1L), symmetric = TRUE)
        loading <- e$vectors
        top <- apply(abs(loading), 2L, which.max)
        loading <- sweep(loading, 2L,
            sign(loading[cbind(top, seq_len(ncol(loading)))]), `*`)
        joint <- joint %*% loading
        eigenvalues <- pmax(e$values, 0)
    }
    dimnames(joint) <- list(cells, paste0("Dim.", seq_len(ncol(joint))))
    block_share <- t(rowsum(loading^2, block_of, reorder = FALSE))
    dimnames(block_share) <- list(colnames(joint), feat_types)

    block_input <- list(
        method = "blockPCA", feat_types = feat_types,
        reduction_methods = reduction_methods,
        reduction_names = reduction_names, dimensions_to_use = dims,
        block_weights = stats::setNames(block_weights, feat_types),
        scale = stats::setNames(scale, feat_types), rotate = rotate,
        block_share = block_share
    )
    for (ft in feat_types) {
        dim_obj <- createDimObj(
            coordinates = joint,
            name = name,
            spat_unit = spat_unit,
            feat_type = ft,
            method = "blockpca",
            reduction = "cells",
            misc = list(eigenvalues = eigenvalues, block_input = block_input)
        )
        gobject <- setDimReduction(
            gobject = gobject,
            x = dim_obj,
            spat_unit = spat_unit,
            feat_type = ft,
            reduction = "cells",
            reduction_method = "blockpca",
            name = name,
            verbose = verbose
        )
    }

    integrated_k <- as.integer(integrated_k)
    if (integrated_k >= length(cells)) {
        stop("[runBlockPCA] integrated_k (", integrated_k, ") must be ",
            "smaller than the number of cells (", length(cells), ")",
            call. = FALSE)
    }
    vmsg(.v = verbose, "runBlockPCA: integrated kNN")
    knn <- .wnn_knn_search(joint, k = integrated_k, engine = engine,
        n_threads = n_threads, dist = TRUE)
    nn_input <- list(
        method = "blockPCA", feat_types = feat_types,
        reduction_methods = "blockpca", reduction_names = name,
        dimensions_to_use = seq_len(ncol(joint)), k = integrated_k,
        engine = engine
    )
    .integrated_nn_store(gobject,
        spat_unit = spat_unit, feat_type = feat_types[1], cells = cells,
        id = knn$id, dist = knn$dist, x = joint,
        network_name = integrated_network_name, snn = integrated_snn,
        snn_prune = snn_prune, snn_name = integrated_snn_name,
        nn_input = nn_input, verbose = verbose
    )
}


.block_weights <- function(block_weights, feat_types) {
    n_mod <- length(feat_types)
    if (is.null(block_weights)) return(rep(1 / n_mod, n_mod))
    checkmate::assert_numeric(block_weights, len = n_mod, finite = TRUE,
        any.missing = FALSE)
    if (any(block_weights <= 0)) {
        stop("[runBlockPCA] `block_weights` must be positive", call. = FALSE)
    }
    if (!is.null(names(block_weights))) {
        if (!setequal(names(block_weights), feat_types)) {
            stop("[runBlockPCA] the names of `block_weights` must be the ",
                "feature types", call. = FALSE)
        }
        block_weights <- block_weights[feat_types]
    }
    unname(block_weights / sum(block_weights))
}
