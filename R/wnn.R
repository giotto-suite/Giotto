#' Multi omics integration with WNN
#'
#' @description
#' Weighted nearest neighbor (WNN) integration of feature types
#' (Hao et al. 2021, *Cell* 184:3573). For each cell and modality, the
#' modality's own neighbors are compared with the other modalities' neighbors
#' at predicting the cell's embedding; a softmax over those scores gives
#' per-cell modality weights, which weight the cell-cell affinities of the
#' integrated neighbors.
#'
#' `graph` chooses how the integrated neighbors are found:
#' * `"dense"` (default): the method of Giotto 4.3.6 and earlier, unchanged.
#'   Two feature types; the weighted affinity of every pair of cells is stored
#'   as a dense N x N matrix (`matrix_result_name`), from which
#'   [runIntegratedUMAP()] builds the integrated kNN. Memory and time grow
#'   with N^2: each N x N matrix takes 8 N^2 bytes (0.75 GB at 10,000 cells,
#'   19 GB at 50,000), and a message says so from 20,000 cells.
#' * `"knn"`: the WNN paper's kernel on k-nearest-neighbor lists only (memory
#'   O(N knn_range), time linear in N), so it scales to whole slides and
#'   GiottoDisk-backed objects. Two or more feature types. The integrated kNN
#'   keeps, for every cell, the `integrated_k` cells with the highest
#'   weighted affinity among the union of each modality's `knn_range` nearest
#'   neighbors, and is stored as a network. Scores are capped and the softmax
#'   is computed in log-sum-exp form, so weights are always finite. Use
#'   `runIntegratedUMAP(graph = "knn")` with it.
#'
#' Arguments marked "`knn` only" are ignored, with a warning, by
#' `graph = "dense"`.
#'
#' @param gobject A Giotto object with a dimension reduction (e.g. PCA) for
#'   each feature type in `feat_types`.
#' @param spat_unit spatial unit
#' @param feat_types feature types to integrate. The integrated kNN is stored
#'   under the first one.
#' @param reduction_methods reduction method of each feature type's embedding.
#'   Default `c("pca", "pca")`.
#' @param reduction_names name of each feature type's embedding. `NULL`
#'   (default): for each feature type, the stored `"<feat_type>.<method>"`
#'   embedding, else `"<method>"`, so both the names of earlier defaults
#'   (`"rna.pca"`) and those `runPCA()` gives (`"pca"`, `"protein.pca"`) are
#'   found. `list_dim_reductions(gobject)` lists the stored names.
#' @param k number of nearest neighbors per modality used for the modality
#'   predictions and the kernel bandwidths. With `graph = "dense"` they are
#'   read from each feature type's stored kNN network.
#' @param integrated_feat_type name under which the results are stored in the
#'   multiomics slot. Default: the feature types joined by `"_"`
#'   (e.g. `"rna_protein"`).
#' @param matrix_result_name name of the stored affinity matrix. `"dense"`:
#'   the N x N matrix of all cell pairs, `NULL` (default) gives
#'   `"theta_weighted_matrix"`. `"knn"`: `NULL` (default) stores none; a name
#'   stores the affinity restricted to the integrated kNN as a sparse N x N
#'   `dgCMatrix` (N x `integrated_k` non-zeros).
#' @param w_names names of the stored weights, one per feature type. `NULL`
#'   (default) gives `"w_<feat_type>"`.
#' @param verbose be verbose
#' @param graph `"dense"` (default, Giotto 4.3.6 and earlier) or `"knn"`
#'   (scales to whole slides); see Description.
#' @param network_names (`knn` only) `NULL` or the name of a stored kNN
#'   network for each feature type whose neighbors are used instead of a new
#'   search (all distances are still computed in the embedding).
#' @param dimensions_to_use (`knn` only) dimensions (indices) of each
#'   embedding to use: one integer vector for all feature types, or a list
#'   with one per feature type. Indices beyond those stored are dropped.
#'   `NULL` (default): the dimensions recorded on the feature type's kNN
#'   network built from the same embedding (`network_names`, else the first
#'   stored kNN network), else all stored dimensions.
#' @param l2_norm (`knn` only) logical. L2-normalize each cell's embedding
#'   before any distance is computed (the WNN paper's default). Default
#'   `TRUE`.
#' @param score_cap (`knn` only) upper limit of the modality scores before the
#'   softmax (default 200, as in the WNN paper). `Inf` leaves them uncapped;
#'   the softmax is stable either way.
#' @param knn_range (`knn` only) number of nearest neighbors (self excluded)
#'   per modality searched as candidates for the integrated kNN.
#' @param integrated_k (`knn` only) number of neighbors (self excluded) in the
#'   integrated kNN. Default `k`.
#' @param integrated_network_name (`knn` only) name of the integrated kNN
#'   network, stored as a `kNN` network under `feat_types[1]`.
#' @param integrated_snn (`knn` only) shared nearest neighbor network of the
#'   integrated kNN to store for clustering, as an `sNN` network named
#'   `integrated_snn_name` under `feat_types[1]`. `"none"` (default) stores no
#'   sNN. `"giotto"` builds it as `createNearestNetwork(type = "sNN")` does:
#'   the kNN edges whose cells share at least 5 neighbors, plus each cell's 3
#'   edges with the most shared neighbors, weighted `1 / (1 + distance)`.
#'   `"jaccard"` links every two cells whose neighbor sets overlap, weighted
#'   by their Jaccard index (shared / (2 `integrated_k` - shared)), and drops
#'   links below `snn_prune`.
#' @param snn_prune (`knn` only) smallest Jaccard index kept with
#'   `integrated_snn = "jaccard"`. Default 1/15.
#' @param integrated_snn_name (`knn` only) name of the sNN network. `NULL`
#'   (default): `integrated_network_name` with `"kNN"` replaced by `"sNN"`.
#' @param engine (`knn` only) nearest neighbor search: `"dbscan"` (exact,
#'   default) or `"hnsw"` (approximate, multithreaded; see
#'   [GiottoClass::hnswKNN()]). The search for `knn_range` neighbors dominates
#'   the run time on large data, and `"hnsw"` is several times faster there:
#'   use it for whole slides.
#' @param n_threads (`knn` only) threads for `engine = "hnsw"`. `NULL` uses
#'   `GiottoUtils::determine_cores()`.
#' @param add_to_metadata (`knn` only) logical. Also add the weights as cell
#'   metadata columns (named by `w_names`) for every feature type in
#'   `feat_types`. Default `FALSE`.
#'
#' @details
#' `graph = "knn"`, for modality m with embedding X_m and k neighbors N_m(i)
#' of cell i: rho_m(i) is the distance to the nearest neighbor; the
#' within-modality prediction is the mean of X_m over N_m(i), the
#' cross-modality prediction from modality o the mean of X_m over N_o(i). The
#' bandwidth b_m(i) is the mean distance (minus rho_m(i)) to the k + 1 cells
#' of lowest non-zero Jaccard similarity in the shared-neighbor graph, as in
#' the WNN paper. theta = exp(-max(d - rho, 0) / b) for the distance d between
#' X_m(i) and each prediction, and S_m,o = theta_within / (theta_cross,o +
#' 1e-4), capped to `[0, score_cap]`. The weight of m is
#' sum_o exp(S_m,o) / sum_m' sum_o exp(S_m',o). A zero bandwidth (identical
#' profiles) gives theta = 1 when the distance equals rho and 0 otherwise.
#' The integrated affinity of a candidate pair is
#' T(a, b) = sum_m w_m(a) exp(-max(D_m(a, b) - rho_m(a), 0) / b_m(a)). The
#' integrated kNN keeps the top `integrated_k` by T (ties by cell order); its
#' edges carry `distance = sqrt(min(max((1 - T) / 2, 0), 1))`,
#' `weight = 1 / (1 + distance)` and `affinity = T`.
#'
#' @returns The Giotto object with, in the multiomics slot under
#'   `integrated_feat_type` and `"WNN"`, one weight vector per feature type
#'   (`w_names`). `"dense"`: also the N x N affinity matrix
#'   (`matrix_result_name`). `"knn"`: the integrated kNN network
#'   (`integrated_network_name`) under `feat_types[1]`, with
#'   `integrated_snn` its sNN (`integrated_snn_name`); and, with
#'   `add_to_metadata = TRUE`, the weights as cell metadata.
#' @seealso [runIntegratedUMAP()]; with `graph = "knn"`, cluster with
#'   `doLeidenCluster(feat_type = feat_types[1], nn_network_to_use = "kNN",
#'   network_name = "integrated_kNN")`, or on the sNN made with
#'   `integrated_snn = "giotto"`.
#' @md
#' @export
runWNN <- function(
        gobject,
        spat_unit = "cell",
        feat_types = c("rna", "protein"),
        reduction_methods = c("pca", "pca"),
        reduction_names = NULL,
        k = 20,
        integrated_feat_type = NULL,
        matrix_result_name = NULL,
        w_names = NULL,
        verbose = FALSE,
        graph = c("dense", "knn"),
        network_names = NULL,
        dimensions_to_use = NULL,
        l2_norm = TRUE,
        score_cap = 200,
        knn_range = 200,
        integrated_k = k,
        integrated_network_name = "integrated_kNN",
        integrated_snn = c("none", "giotto", "jaccard"),
        snn_prune = 1 / 15,
        integrated_snn_name = NULL,
        engine = c("dbscan", "hnsw"),
        n_threads = NULL,
        add_to_metadata = FALSE) {
    checkmate::assert_class(gobject, "giotto")
    call_args <- names(match.call())
    graph <- match.arg(graph)
    checkmate::assert_string(spat_unit)
    checkmate::assert_character(feat_types, min.len = 2L, unique = TRUE,
        any.missing = FALSE)
    n_mod <- length(feat_types)
    if (is.null(w_names)) w_names <- paste0("w_", feat_types)
    checkmate::assert_character(w_names, len = n_mod, unique = TRUE,
        any.missing = FALSE)
    reduction_methods <- .wnn_per_mod(reduction_methods, feat_types,
        "reduction_methods")
    .wnn_check_reductions(gobject, spat_unit, feat_types, "runWNN")
    if (is.null(reduction_names)) {
        reduction_names <- .wnn_reduction_names(gobject, spat_unit,
            feat_types, reduction_methods)
    }
    reduction_names <- .wnn_per_mod(reduction_names, feat_types,
        "reduction_names")

    if (identical(graph, "dense")) {
        .wnn_ignored(call_args, .wnn_knn_args, "runWNN")
        if (n_mod != 2L) {
            stop("[runWNN] graph = \"dense\" integrates two feature types; ",
                "use graph = \"knn\" for more", call. = FALSE)
        }
        .wnn_dense_notice(length(spatIDs(gobject, spat_unit = spat_unit)))
        return(.wnn_dense(gobject,
            spat_unit = spat_unit, feat_types = feat_types,
            reduction_methods = reduction_methods,
            reduction_names = reduction_names, k = k,
            integrated_feat_type = integrated_feat_type,
            matrix_result_name = matrix_result_name, w_names = w_names,
            verbose = verbose
        ))
    }

    engine <- match.arg(engine)
    integrated_snn <- match.arg(integrated_snn)
    if (!is.null(network_names)) {
        network_names <- .wnn_per_mod(network_names, feat_types,
            "network_names")
    }
    checkmate::assert_flag(l2_norm)
    checkmate::assert_count(k, positive = TRUE)
    checkmate::assert_count(knn_range, positive = TRUE)
    checkmate::assert_count(integrated_k, positive = TRUE)
    checkmate::assert_number(score_cap, lower = 0)
    checkmate::assert_string(integrated_network_name)
    checkmate::assert_number(snn_prune, lower = 0, upper = 1)
    if (is.null(integrated_snn_name)) {
        integrated_snn_name <- sub("kNN", "sNN", integrated_network_name,
            fixed = TRUE)
    }
    checkmate::assert_string(integrated_snn_name)
    checkmate::assert_flag(add_to_metadata)
    k <- as.integer(k)
    if (k < 2L) stop("[runWNN] k must be at least 2", call. = FALSE)
    if (integrated_k > knn_range) {
        stop("[runWNN] integrated_k (", integrated_k, ") cannot exceed ",
            "knn_range (", knn_range, ")", call. = FALSE)
    }
    if (is.null(integrated_feat_type)) {
        integrated_feat_type <- paste(feat_types, collapse = "_")
    }

    inp <- .wnn_inputs(gobject,
        spat_unit = spat_unit, feat_types = feat_types,
        reduction_methods = reduction_methods,
        reduction_names = reduction_names, network_names = network_names,
        dimensions_to_use = dimensions_to_use, l2_norm = l2_norm,
        verbose = verbose
    )
    cells <- inp$cells
    n_cells <- length(cells)
    if (n_cells <= k) {
        stop("[runWNN] k (", k, ") must be smaller than the number of ",
            "cells (", n_cells, ")", call. = FALSE)
    }
    knn_range <- min(as.integer(knn_range), n_cells - 1L)
    integrated_k <- min(as.integer(integrated_k), knn_range)
    vmsg(.v = verbose, "runWNN: searching neighbors")
    nb <- cand <- list()
    for (ft in feat_types) {
        need_search <- is.null(inp$nb[[ft]])
        kk <- max(knn_range, if (need_search) k else 0L)
        knn <- .wnn_knn_search(inp$X[[ft]], k = kk, engine = engine,
            n_threads = n_threads)
        cand[[ft]] <- knn$id[, seq_len(knn_range), drop = FALSE]
        nb[[ft]] <- if (need_search) {
            knn$id[, seq_len(k), drop = FALSE]
        } else {
            .wnn_closest(inp$X[[ft]], inp$nb[[ft]], k)
        }
        rm(knn)
    }

    vmsg(.v = verbose, "runWNN: modality weights")
    mw <- .wnn_modality_weights(inp$X, nb, score_cap = score_cap)

    vmsg(.v = verbose, "runWNN: integrated kNN")
    ik <- .wnn_integrated_knn(inp$X, cand,
        w = mw$w, bw = mw$bw, rho = mw$rho, k = integrated_k
    )
    rownames(ik$id) <- rownames(ik$dist) <- cells

    for (i in seq_len(n_mod)) {
        w_i <- stats::setNames(mw$w[[i]], cells)
        gobject <- setMultiomics(
            gobject = gobject,
            x = stats::setNames(list(w_i), feat_types[i]),
            spat_unit = spat_unit,
            feat_type = integrated_feat_type,
            integration_method = "WNN",
            result_name = w_names[i],
            verbose = verbose
        )
    }
    if (!is.null(matrix_result_name)) {
        tw <- Matrix::sparseMatrix(
            i = rep(seq_len(n_cells), integrated_k),
            j = as.vector(ik$id), x = as.vector(ik$affinity),
            dims = c(n_cells, n_cells), dimnames = list(cells, cells)
        )
        gobject <- setMultiomics(gobject,
            x = tw, spat_unit = spat_unit, feat_type = integrated_feat_type,
            integration_method = "WNN", result_name = matrix_result_name,
            verbose = verbose
        )
    }

    nn_input <- list(
        method = "WNN", feat_types = feat_types,
        reduction_methods = reduction_methods,
        reduction_names = reduction_names,
        dimensions_to_use = inp$dims, network_names = inp$network_names,
        l2_norm = l2_norm, score_cap = score_cap,
        k = integrated_k, k_modality = k, knn_range = knn_range,
        engine = engine, integrated_feat_type = integrated_feat_type
    )
    gobject <- .integrated_nn_store(gobject,
        spat_unit = spat_unit, feat_type = feat_types[1], cells = cells,
        id = ik$id, dist = ik$dist, affinity = ik$affinity, x = inp$X[[1L]],
        network_name = integrated_network_name, snn = integrated_snn,
        snn_prune = snn_prune, snn_name = integrated_snn_name,
        nn_input = nn_input, verbose = verbose
    )

    if (add_to_metadata) {
        meta <- data.table::data.table(cell_ID = cells)
        for (i in seq_len(n_mod)) {
            data.table::set(meta, j = w_names[i], value = mw$w[[i]])
        }
        for (ft in feat_types) {
            gobject <- addCellMetadata(gobject,
                spat_unit = spat_unit, feat_type = ft,
                new_metadata = meta, by_column = TRUE,
                column_cell_ID = "cell_ID"
            )
        }
    }

    gobject
}


#' Run integrated UMAP
#'
#' @description
#' UMAP of the integrated neighbors made by [runWNN()], with the same `graph`:
#' * `"dense"` (default, Giotto 4.3.6 and earlier, unchanged): builds the
#'   integrated kNN from the N x N affinity matrix of `runWNN(graph =
#'   "dense")` unless one is stored (`force = TRUE` rebuilds it), and embeds
#'   it.
#' * `"knn"`: uwot embeds the integrated kNN network of `runWNN(graph =
#'   "knn")` or [runBlockPCA()] directly (`nn_method`), so no data matrix is
#'   needed.
#'
#' Arguments marked "`dense` only" or "`knn` only" are ignored, with a
#' warning, by the other `graph`.
#'
#' @param gobject A giotto object
#' @param spat_unit spatial unit
#' @param feat_types feature types that were integrated. The UMAP is stored
#'   under each of them; the integrated kNN is read from the first.
#' @param k number of neighbors (self excluded). With `graph = "knn"`, taken
#'   from the integrated kNN (at most `runWNN(integrated_k =)`).
#' @param spread UMAP param: spread
#' @param min_dist UMAP param: min_dist
#' @param integrated_feat_type (`dense` only) integrated feature type
#'   (e.g. 'rna_protein')
#' @param integration_method (`dense` only) multiomics integration method
#'   used. Default = 'WNN'
#' @param matrix_result_name (`dense` only) Default = 'theta_weighted_matrix'
#' @param force (`dense` only) force calculation of integrated kNN.
#'   Default = FALSE
#' @param seed seed. Default 1234
#' @param graph `"dense"` (default) or `"knn"`, as used in [runWNN()].
#' @param network_name (`knn` only) name of the integrated kNN network made
#'   by [runWNN()] or [runBlockPCA()]. Default `"integrated_kNN"`.
#' @param name (`knn` only) name of the stored UMAP. Default
#'   `"integrated.umap"` (the name `graph = "dense"` uses).
#' @param ... additional parameters for [uwot::umap()] (e.g. `n_threads`)
#'
#' @returns A Giotto object with the integrated UMAP stored as a `umap`
#'   reduction for every feature type in `feat_types`.
#' @md
#' @export
runIntegratedUMAP <- function(
        gobject,
        spat_unit = "cell",
        feat_types = c("rna", "protein"),
        integrated_feat_type = NULL,
        integration_method = "WNN",
        matrix_result_name = "theta_weighted_matrix",
        k = 20,
        spread = 5,
        min_dist = 0.01,
        force = FALSE,
        seed = 1234,
        graph = c("dense", "knn"),
        network_name = "integrated_kNN",
        name = "integrated.umap",
        ...) {
    checkmate::assert_class(gobject, "giotto")
    graph <- match.arg(graph)
    call_args <- names(match.call())
    if (identical(graph, "dense")) {
        .wnn_ignored(call_args, c("network_name", "name"),
            "runIntegratedUMAP")
        return(.wnn_umap_dense(gobject,
            spat_unit = spat_unit, feat_types = feat_types,
            integrated_feat_type = integrated_feat_type,
            integration_method = integration_method,
            matrix_result_name = matrix_result_name, k = k, spread = spread,
            min_dist = min_dist, force = force, seed = seed, ...
        ))
    }
    .wnn_ignored(call_args, c("integrated_feat_type", "integration_method",
        "matrix_result_name", "force"), "runIntegratedUMAP", "knn")
    .wnn_check_reductions(gobject, spat_unit, feat_types,
        "runIntegratedUMAP")
    checkmate::assert_count(k, positive = TRUE)
    checkmate::assert_string(network_name)
    checkmate::assert_string(name)
    nn <- .wnn_stored_integrated_knn(gobject,
        spat_unit = spat_unit, feat_type = feat_types[1],
        network_name = network_name
    )
    if (k > ncol(nn$id)) {
        stop("[runIntegratedUMAP] k (", k, ") is larger than the ",
            "integrated kNN (", ncol(nn$id), " neighbors per cell); rerun ",
            "runWNN() or runBlockPCA() with a larger integrated_k",
            call. = FALSE)
    }
    id <- nn$id[, seq_len(k), drop = FALSE]
    dist <- nn$dist[, seq_len(k), drop = FALSE]
    cells <- rownames(nn$id)

    GiottoUtils::local_seed(seed)
    integrated_umap <- uwot::umap(
        X = NULL,
        n_neighbors = k,
        nn_method = GiottoClass::nnToUwot(list(id = id, dist = dist)),
        spread = spread,
        min_dist = min_dist,
        ...
    )
    dimnames(integrated_umap) <- list(cells,
        paste0("Dim.", seq_len(ncol(integrated_umap))))

    for (ft in feat_types) {
        umap_object <- createDimObj(
            coordinates = integrated_umap,
            spat_unit = spat_unit,
            feat_type = ft,
            method = "umap",
            name = name
        )
        gobject <- setDimReduction(
            gobject = gobject,
            x = umap_object,
            spat_unit = spat_unit,
            feat_type = ft,
            reduction_method = "umap",
            name = name
        )
    }

    gobject
}


# arguments of runWNN() that only graph = "knn" uses
.wnn_knn_args <- c("network_names", "dimensions_to_use", "l2_norm",
    "score_cap", "knn_range", "integrated_k", "integrated_network_name",
    "integrated_snn", "snn_prune", "integrated_snn_name", "engine",
    "n_threads", "add_to_metadata")

.wnn_ignored <- function(call_args, args, fun, graph = "dense") {
    set <- intersect(call_args, args)
    if (length(set)) {
        warning(GiottoUtils::wrap_txt(sprintf(
            "[%s] ignored with graph = \"%s\": %s", fun, graph,
            paste(set, collapse = ", "))), call. = FALSE)
    }
    invisible(set)
}

.wnn_dense_notice <- function(n_cells, min_cells = 20000L) {
    if (n_cells < min_cells || !getOption("giotto.warn_wnn_dense", TRUE)) {
        return(invisible(FALSE))
    }
    vmsg(.v = TRUE, sprintf(paste(
        "runWNN(graph = \"dense\") builds N x N matrices of all cell pairs:",
        "%.1f GB each for %s cells. For large data use graph = \"knn\" in",
        "runWNN() and runIntegratedUMAP()."), 8 * n_cells^2 / 2^30,
        format(n_cells, big.mark = ",")))
    invisible(TRUE)
}

.wnn_reduction_names <- function(gobject, spat_unit, feat_types,
    reduction_methods) {
    dr <- methods::slot(gobject, "dimension_reduction")$cells[[spat_unit]]
    vapply(seq_along(feat_types), function(i) {
        ft <- feat_types[i]
        rm <- reduction_methods[i]
        cand <- c(paste0(ft, ".", rm), rm)
        hit <- intersect(cand, names(dr[[ft]][[rm]]))
        if (length(hit)) hit[1L] else cand[if (ft == "rna") 2L else 1L]
    }, character(1L))
}

.wnn_check_reductions <- function(gobject, spat_unit, feat_types, fun) {
    have <- names(methods::slot(gobject, "dimension_reduction")$cells[[
        spat_unit]])
    absent <- setdiff(feat_types, have)
    if (length(absent)) {
        stop(GiottoUtils::wrap_txt(paste0("[", fun, "]"),
            "no dimension reduction for feature type(s):",
            paste(absent, collapse = ", ")), call. = FALSE)
    }
}

.wnn_per_mod <- function(x, feat_types, what, fun = "runWNN") {
    if (length(unique(x)) == 1L) x <- rep(x[1L], length(feat_types))
    if (length(x) != length(feat_types)) {
        stop("[", fun, "] `", what, "` needs one value per feature type (",
            length(feat_types), ")", call. = FALSE)
    }
    as.character(x)
}

.wnn_per_mod_dims <- function(dimensions_to_use, feat_types,
    fun = "runWNN") {
    n_mod <- length(feat_types)
    if (is.null(dimensions_to_use)) return(vector("list", n_mod))
    if (!is.list(dimensions_to_use)) {
        return(rep(list(dimensions_to_use), n_mod))
    }
    if (length(dimensions_to_use) != n_mod) {
        stop("[", fun, "] `dimensions_to_use` needs one entry per feature ",
            "type (", n_mod, ")", call. = FALSE)
    }
    if (!is.null(names(dimensions_to_use))) {
        if (!setequal(names(dimensions_to_use), feat_types)) {
            stop("[", fun, "] the names of `dimensions_to_use` must be the ",
                "feature types", call. = FALSE)
        }
        dimensions_to_use <- dimensions_to_use[feat_types]
    }
    dimensions_to_use
}

.wnn_inputs <- function(gobject, spat_unit, feat_types, reduction_methods,
    reduction_names, network_names, dimensions_to_use, l2_norm,
    verbose = NULL) {
    n_mod <- length(feat_types)
    if (!is.null(dimensions_to_use)) {
        dimensions_to_use <- .wnn_per_mod_dims(dimensions_to_use, feat_types)
    }

    X_native <- X <- nb <- dims <- list()
    net_used <- character(n_mod)
    cells <- NULL
    for (i in seq_len(n_mod)) {
        ft <- feat_types[i]
        emb <- getDimReduction(gobject,
            spat_unit = spat_unit, feat_type = ft, reduction = "cells",
            reduction_method = reduction_methods[i],
            name = reduction_names[i], output = "matrix"
        )
        if (is.null(rownames(emb)) || anyDuplicated(rownames(emb))) {
            stop("[runWNN] the embedding of ", ft, " needs unique cell ",
                "IDs as row names", call. = FALSE)
        }

        net_name <- if (!is.null(network_names)) network_names[i]
        d <- if (!is.null(dimensions_to_use)) {
            dimensions_to_use[[i]]
        } else {
            .wnn_recorded_dims(gobject, spat_unit, ft, net_name,
                reduction_methods[i], reduction_names[i])
        }
        all_d <- seq_len(ncol(emb))
        d <- if (is.null(d)) all_d else intersect(as.integer(d), all_d)
        if (!length(d)) {
            stop("[runWNN] no stored dimensions of ", ft, " are selected",
                call. = FALSE)
        }
        vmsg(.v = verbose, sprintf("runWNN: %s uses dimensions %s", ft,
            .wnn_range_txt(d)))
        emb <- emb[, d, drop = FALSE]
        storage.mode(emb) <- "double"
        bad <- rowSums(!is.finite(emb)) > 0L
        if (any(bad)) {
            stop(sprintf(paste0("[runWNN] the embedding of %s has NA or ",
                "infinite values for %d cells"), ft, sum(bad)), call. = FALSE)
        }
        if (l2_norm && length(d) < 2L) {
            stop(sprintf(paste0("[runWNN] l2_norm needs at least 2 ",
                "dimensions of %s (got %d); use more dimensions or ",
                "l2_norm = FALSE"), ft, length(d)), call. = FALSE)
        }
        if (l2_norm) emb <- .wnn_l2(emb)
        X_native[[ft]] <- emb

        if (is.null(cells)) {
            cells <- rownames(emb)
        } else if (!setequal(cells, rownames(emb))) {
            n_miss <- sum(!cells %in% rownames(emb))
            n_extra <- sum(!rownames(emb) %in% cells)
            stop(sprintf(paste0("[runWNN] %s and %s embed different cells ",
                "(%d missing from %s, %d only in %s); subset them to the ",
                "same cells first"), feat_types[1], ft, n_miss, ft, n_extra,
                ft), call. = FALSE)
        }
        X[[ft]] <- if (identical(rownames(emb), cells)) {
            emb
        } else {
            emb[cells, , drop = FALSE]
        }
        dims[[ft]] <- d
        if (!is.null(net_name)) {
            nb[[ft]] <- .wnn_network_ids(gobject, spat_unit, ft, net_name,
                cells)
            net_used[i] <- net_name
        }
    }
    list(X = X, X_native = X_native, nb = nb, cells = cells, dims = dims,
        network_names = net_used)
}

.wnn_range_txt <- function(d) {
    if (length(d) > 1L && all(diff(d) == 1L)) {
        paste0(d[1], ":", d[length(d)])
    } else {
        paste(d, collapse = ",")
    }
}

.wnn_l2 <- function(x) {
    nrm <- sqrt(rowSums(x^2))
    out <- x / nrm
    out[!is.finite(out)] <- 0
    out
}

.wnn_knn_nets <- function(gobject, spat_unit, feat_type) {
    nets <- methods::slot(gobject, "nn_network")[[spat_unit]][[feat_type]][[
        "kNN"]]
    is_wnn <- vapply(nets, function(x) {
        identical(.wnn_nn_input(x)$method, "WNN")
    }, logical(1L))
    nets[!is_wnn]
}

.wnn_nn_input <- function(x) {
    tryCatch(x@misc$nn_input, error = function(e) NULL)
}

.wnn_recorded_dims <- function(gobject, spat_unit, feat_type, net_name,
    reduction_method, reduction_name) {
    nets <- .wnn_knn_nets(gobject, spat_unit, feat_type)
    if (!is.null(net_name)) nets <- nets[names(nets) == net_name]
    for (obj in nets) {
        rec <- .wnn_nn_input(obj)
        if (!is.null(rec$dimensions_to_use) &&
            identical(rec$dim_reduction_to_use, reduction_method) &&
            identical(rec$dim_reduction_name, reduction_name)) {
            return(rec$dimensions_to_use)
        }
    }
    NULL
}

.wnn_network_ids <- function(gobject, spat_unit, feat_type, name, cells) {
    dt <- getNearestNetwork(gobject,
        spat_unit = spat_unit, feat_type = feat_type, nn_type = "kNN",
        name = name, output = "data.table"
    )
    dt <- data.table::as.data.table(dt)
    from_col <- intersect(c("from", "from_id"), names(dt))[1L]
    to_col <- intersect(c("to", "to_id"), names(dt))[1L]
    f <- match(dt[[from_col]], cells)
    t <- match(dt[[to_col]], cells)
    if (anyNA(f) || anyNA(t)) {
        stop(sprintf(paste0("[runWNN] the kNN network '%s' of %s covers ",
            "cells that are not in the embedding"), name, feat_type),
            call. = FALSE)
    }
    cnt <- tabulate(f, nbins = length(cells))
    if (any(cnt != cnt[1L]) || cnt[1L] == 0L) {
        stop(sprintf(paste0("[runWNN] the kNN network '%s' of %s does not ",
            "give every cell the same number of neighbors (%d-%d)"),
            name, feat_type, min(cnt), max(cnt)), call. = FALSE)
    }
    o <- order(f)
    matrix(t[o], ncol = cnt[1L], byrow = TRUE)
}

.wnn_closest <- function(X, nb, k) {
    if (ncol(nb) < k) {
        stop(sprintf(paste0("[runWNN] the stored kNN network has %d ",
            "neighbors per cell, k = %d needs at least that many"),
            ncol(nb), k), call. = FALSE)
    }
    if (ncol(nb) == k) return(nb)
    dk <- .wnn_nb_dist(X, nb)
    o <- .wnn_row_order(dk, nb)
    o[, seq_len(k), drop = FALSE]
}

.wnn_row_order <- function(val, id) {
    n <- nrow(id)
    ord <- order(rep(seq_len(n), ncol(id)), as.vector(val), as.vector(id),
        method = "radix")
    matrix(as.vector(id)[ord], nrow = n, byrow = TRUE)
}

.wnn_knn_search <- function(x, k, engine = c("dbscan", "hnsw"),
    n_threads = NULL, dist = FALSE) {
    engine <- match.arg(engine)
    res <- if (identical(engine, "dbscan")) {
        dbscan::kNN(x = x, k = k, sort = TRUE)
    } else {
        GiottoClass::hnswKNN(x = x, k = k, n_threads = n_threads)
    }
    id <- res$id
    storage.mode(id) <- "integer"
    dimnames(id) <- NULL
    if (!dist) return(list(id = id))
    d <- res$dist
    dimnames(d) <- NULL
    list(id = id, dist = d)
}

# integrated kNN and, unless `snn` is "none", its sNN; used by runWNN()
# and runBlockPCA()
.integrated_nn_store <- function(gobject, spat_unit, feat_type, cells, id,
    dist, affinity = NULL, x, network_name, snn, snn_prune, snn_name,
    nn_input, verbose = FALSE) {
    k <- ncol(id)
    n_cells <- length(cells)
    edges <- data.table::data.table(
        from = cells[rep(seq_len(n_cells), each = k)],
        to = cells[as.vector(t(id))],
        weight = 1 / (1 + as.vector(t(dist))),
        distance = as.vector(t(dist))
    )
    if (!is.null(affinity)) {
        data.table::set(edges, j = "affinity", value = as.vector(t(affinity)))
    }
    net <- igraph::graph_from_data_frame(edges,
        directed = TRUE,
        vertices = data.frame(name = cells, stringsAsFactors = FALSE)
    )
    nn_obj <- createNearestNetObj(
        name = network_name, network = net, nn_type = "kNN",
        spat_unit = spat_unit, feat_type = feat_type,
        misc = list(nn_input = nn_input)
    )
    gobject <- setNearestNetwork(gobject, nn_obj,
        spat_unit = spat_unit, feat_type = feat_type, nn_type = "kNN",
        name = network_name, verbose = verbose
    )
    if (identical(snn, "none")) return(gobject)
    snn_input <- c(nn_input, list(type = "sNN", snn = snn))
    if (identical(snn, "jaccard")) {
        snn_net <- .wnn_jaccard_snn(id, cells, prune = snn_prune)
        snn_input$prune <- snn_prune
    } else {
        knn <- structure(list(dist = unname(dist), id = unname(id), k = k,
            sort = TRUE), class = c("kNN", "NN"))
        snn_net <- GiottoClass::createNetwork(x,
            GiottoClass::sNNNetworkParam(k = k, output = "igraph"),
            node_ids = cells, nn_network = knn, verbose = FALSE
        )
        snn_input$minimum_shared <- 5L
        snn_input$top_shared <- 3L
    }
    snn_obj <- createNearestNetObj(
        name = snn_name, network = snn_net, nn_type = "sNN",
        spat_unit = spat_unit, feat_type = feat_type,
        misc = list(nn_input = snn_input)
    )
    setNearestNetwork(gobject, snn_obj,
        spat_unit = spat_unit, feat_type = feat_type, nn_type = "sNN",
        name = snn_name, verbose = verbose
    )
}

.wnn_nb_dist <- function(X, nb) {
    out <- matrix(0, nrow(nb), ncol(nb))
    for (j in seq_len(ncol(nb))) {
        out[, j] <- sqrt(rowSums((X - X[nb[, j], , drop = FALSE])^2))
    }
    out
}

.wnn_pred_dist <- function(X, nb) {
    acc <- X[nb[, 1L], , drop = FALSE]
    if (ncol(nb) > 1L) {
        for (j in 2:ncol(nb)) acc <- acc + X[nb[, j], , drop = FALSE]
    }
    sqrt(rowSums((X - acc / ncol(nb))^2))
}

.wnn_pair_dist <- function(X, i, j, batch = 131072L) {
    n <- length(i)
    out <- numeric(n)
    if (!n) return(out)
    for (s in seq(1L, n, by = batch)) {
        r <- s:min(s + batch - 1L, n)
        out[r] <- sqrt(rowSums((X[i[r], , drop = FALSE] -
            X[j[r], , drop = FALSE])^2))
    }
    out
}

.wnn_relu <- function(x) {
    (x + abs(x)) * 0.5
}

.wnn_kernel <- function(num, b) {
    out <- exp(-num / b)
    zero <- !(b > 0)
    if (any(zero)) {
        z <- if (length(b) == 1L) rep(TRUE, length(num)) else
            rep_len(zero, length(num))
        out[z] <- as.numeric(num[z] == 0)
    }
    out
}

.wnn_row_min <- function(M) {
    do.call(pmin, lapply(seq_len(ncol(M)), function(j) M[, j]))
}

.wnn_modality_weights <- function(X, nb, score_cap = 200, epsilon = 1e-4) {
    mods <- names(X)
    n <- nrow(X[[1L]])
    rho <- bw <- S <- list()
    for (m in mods) {
        dk <- .wnn_nb_dist(X[[m]], nb[[m]])
        rho_m <- .wnn_row_min(dk)
        d_w <- .wnn_pred_dist(X[[m]], nb[[m]])
        d_c <- vapply(setdiff(mods, m), function(o) {
            .wnn_pred_dist(X[[m]], nb[[o]])
        }, numeric(n))
        if (!is.matrix(d_c)) d_c <- matrix(d_c, nrow = n)
        bw_m <- .wnn_sigma_snn(X[[m]], nb[[m]], rho_m)
        th_w <- .wnn_kernel(.wnn_relu(d_w - rho_m), bw_m)
        th_c <- .wnn_kernel(.wnn_relu(d_c - rho_m), bw_m)
        S[[m]] <- th_w / (th_c + epsilon)
        rho[[m]] <- rho_m
        bw[[m]] <- bw_m
    }
    Sall <- do.call(cbind, S)
    Sall <- pmin(pmax(Sall, 0), score_cap)
    shift <- do.call(pmax, lapply(seq_len(ncol(Sall)), function(j) Sall[, j]))
    E <- exp(Sall - shift)
    owner <- rep(seq_along(mods), vapply(S, NCOL, integer(1L)))
    num <- vapply(seq_along(mods), function(i) {
        rowSums(E[, owner == i, drop = FALSE])
    }, numeric(n))
    if (!is.matrix(num)) num <- matrix(num, nrow = n)
    W <- num / rowSums(num)
    w <- lapply(seq_along(mods), function(i) W[, i])
    names(w) <- mods
    list(w = w, rho = rho, bw = bw, score = S)
}

.wnn_sigma_snn <- function(X, nb, rho, n_far = ncol(nb) + 1L,
    chunk = NULL) {
    n <- nrow(X)
    k_snn <- ncol(nb) + 1L
    At <- Matrix::sparseMatrix(
        i = c(seq_len(n), as.vector(nb)), j = rep(seq_len(n), k_snn),
        x = 1, dims = c(n, n)
    )
    A <- Matrix::t(At)
    chunk <- chunk %null% max(256L, floor(4e6 / k_snn^2))
    out <- numeric(n)
    for (s in seq(1L, n, by = chunk)) {
        rows <- s:min(s + chunk - 1L, n)
        nr <- length(rows)
        sh <- A %*% At[, rows, drop = FALSE]
        n_col <- diff(sh@p)
        col <- rep.int(seq_len(nr), n_col)
        shared <- as.integer(round(sh@x))
        tab <- matrix(tabulate((col - 1L) * k_snn + shared,
            nbins = nr * k_snn), nrow = k_snn)
        cut <- colSums(MatrixGenerics::colCumsums(tab) <
            pmin(n_far, n_col)) + 1L
        keep <- shared <= cut[col]
        col <- col[keep]
        d <- .wnn_relu(.wnn_pair_dist(X, rows[col], sh@i[keep] + 1L) -
            rho[rows[col]])
        rm(sh, keep, shared)
        o <- order(col, -d, method = "radix")
        col <- col[o]
        d <- d[o]
        top <- seq_along(col) - cumsum(c(1L, tabulate(col, nbins = nr)))[
            col] < n_far
        out[rows] <- rowsum(d[top], col[top], reorder = TRUE)[, 1L] /
            tabulate(col[top], nbins = nr)
    }
    out
}

.wnn_integrated_knn <- function(X, cand, w, bw, rho, k, chunk = NULL) {
    mods <- names(X)
    n <- nrow(X[[1L]])
    width <- sum(vapply(cand, ncol, integer(1L)))
    chunk <- chunk %null% max(64L, floor(2e6 / width))
    id <- matrix(NA_integer_, n, k)
    affinity <- matrix(NA_real_, n, k)
    for (s in seq(1L, n, by = chunk)) {
        rows <- s:min(s + chunk - 1L, n)
        nr <- length(rows)
        ci <- rep.int(seq_len(nr), width)
        j <- unlist(lapply(cand, function(cm) {
            as.vector(cm[rows, , drop = FALSE])
        }), use.names = FALSE)
        o <- order(ci, j, method = "radix")
        ci <- ci[o]
        j <- j[o]
        keep <- c(TRUE, ci[-1L] != ci[-length(ci)] | j[-1L] != j[-length(j)])
        keep <- keep & j != rows[ci]
        ci <- ci[keep]
        j <- j[keep]
        i <- rows[ci]
        tt <- numeric(length(j))
        for (m in mods) {
            dd <- .wnn_pair_dist(X[[m]], i, j)
            kern <- .wnn_kernel(.wnn_relu(dd - rho[[m]][i]), bw[[m]][i])
            tt <- tt + w[[m]][i] * kern
        }
        o <- order(ci, -tt, j, method = "radix")
        ci <- ci[o]
        top <- seq_along(ci) - cumsum(c(1L, tabulate(ci, nbins = nr)))[ci] < k
        if (sum(top) != nr * k) {
            stop("[runWNN] fewer than ", k, " candidate neighbors for some ",
                "cells; increase knn_range", call. = FALSE)
        }
        id[rows, ] <- matrix(j[o][top], ncol = k, byrow = TRUE)
        affinity[rows, ] <- matrix(tt[o][top], ncol = k, byrow = TRUE)
    }
    dist <- sqrt(pmin(pmax((1 - affinity) / 2, 0), 1))
    list(id = id, dist = dist, affinity = affinity)
}

.wnn_jaccard_snn <- function(id, cells, prune = 1 / 15, chunk = NULL) {
    n <- nrow(id)
    k <- ncol(id)
    A <- Matrix::sparseMatrix(i = rep(seq_len(n), k), j = as.vector(id),
        x = 1, dims = c(n, n))
    At <- Matrix::t(A)
    chunk <- chunk %null% max(256L, floor(4e6 / k^2))
    edges <- list()
    for (s in seq(1L, n, by = chunk)) {
        rows <- s:min(s + chunk - 1L, n)
        sh <- A %*% At[, rows, drop = FALSE]
        to <- rep.int(rows, diff(sh@p))
        from <- sh@i + 1L
        shared <- sh@x
        jac <- shared / (2 * k - shared)
        keep <- from < to & jac >= prune
        edges[[length(edges) + 1L]] <- data.table::data.table(
            from = from[keep], to = to[keep], shared = shared[keep],
            jac = jac[keep])
    }
    edges <- data.table::rbindlist(edges)
    igraph::graph_from_data_frame(
        data.frame(from = cells[edges$from], to = cells[edges$to],
            weight = edges$jac, distance = 1 - edges$jac,
            shared = edges$shared),
        directed = FALSE,
        vertices = data.frame(name = cells, stringsAsFactors = FALSE)
    )
}

.wnn_stored_integrated_knn <- function(gobject, spat_unit, feat_type,
    network_name) {
    cells <- spatIDs(gobject, spat_unit = spat_unit)
    nets <- methods::slot(gobject, "nn_network")[[spat_unit]][[feat_type]][[
        "kNN"]]
    if (is.null(nets[[network_name]])) {
        stop(sprintf(paste0("[runIntegratedUMAP] no integrated kNN '%s' ",
            "under %s; run runWNN() or runBlockPCA() first"), network_name,
            feat_type),
            call. = FALSE)
    }
    nn <- tryCatch(.umap_stored_knn(gobject,
        spat_unit = spat_unit, feat_type = feat_type, cell_ids = cells,
        n_neighbors = NULL, nn_type = "kNN", network_name = network_name,
        required = TRUE
    ), error = function(e) {
        stop(sprintf(paste0("[runIntegratedUMAP] the integrated kNN '%s' ",
            "does not match the object's cells; rerun runWNN() or ",
            "runBlockPCA() after subsetting (%s)"), network_name,
            conditionMessage(e)),
            call. = FALSE)
    })
    list(id = nn$id, dist = nn$dist)
}


# runWNN() of Giotto 4.3.6 (graph = "dense"), unchanged
.wnn_dense <- function(
        gobject,
        spat_unit = "cell",
        feat_types = c("rna", "protein"),
        reduction_methods = c("pca", "pca"),
        reduction_names = c("rna.pca", "protein.pca"),
        k = 20,
        integrated_feat_type = NULL,
        matrix_result_name = NULL,
        w_names = c(NULL, NULL),
        verbose = FALSE) {
    # validate Giotto object
    if (!inherits(gobject, "giotto")) {
        stop("gobject needs to be a giotto object")
    }

    # validate feat_types
    for (feat_type in feat_types) {
        if (!feat_type %in% names(
            slot(gobject, "dimension_reduction")$cells[[spat_unit]]
        )) {
            stop(paste(feat_type, " and their dimension reductions must exist in the Giotto object"))
        }
    }

    # extract kNN and PCA
    kNN_list <- list()
    pca_list <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        kNN_list[[feat_type]] <- getNearestNetwork(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            nn_type = "kNN",
            output = "igraph"
        )

        pca_list[[feat_type]] <- getDimReduction(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction = "cells",
            reduction_method = reduction_methods[i],
            name = reduction_names[i],
            output = "matrix"
        )
    }


    ## get cell names
    cell_names <- spatIDs(gobject, spat_unit = spat_unit)

    ######################## distances calculation ############################

    cell_distances <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        cell_distances[[feat_type]] <- list()

        message("Calculating ", feat_type, "-", feat_type, " distance")

        cell_distances[[feat_type]][[feat_type]] <- list()

        for (cell_a in cell_names) {
            my_kNN <- kNN_list[[feat_type]][[cell_a]][[cell_a]]

            cell_distances[[feat_type]][[feat_type]][[cell_a]] <- rep(0, k)
            names(cell_distances[[feat_type]][[feat_type]][[cell_a]]) <- names(my_kNN)

            for (cell_i in names(my_kNN)) {
                dimensions_cell_a_i <- pca_list[[feat_type]][c(cell_a, cell_i), ]
                cell_distances[[feat_type]][[feat_type]][[cell_a]][cell_i] <- sqrt(sum((
                    dimensions_cell_a_i[1, ] - dimensions_cell_a_i[2, ])^2))
            }
        }
    }

    ########################### all cell-cell distances ########################

    all_cell_distances <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        all_cell_distances[[feat_type]] <- list()

        if (verbose) {
            message(paste(
                "Calculating low dimensional cell-cell distances for",
                feat_type
            ))
        }

        pca_distances <- dist(pca_list[[feat_type]])
        all_cell_distances[[feat_type]][[feat_type]] <- as.matrix(pca_distances)
    }

    ######################## within-modality prediction ########################

    if (verbose) message("Calculating within-modality prediction")

    ### predicted feat_type1 feat_type1
    predicted_values <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        predicted_values[[feat_type]] <- list()
        predicted_values[[feat_type]][[feat_type]] <- list()

        for (cell_a in cell_names) {
            dimensions_cell_a <- pca_list[[feat_type]][kNN_list[[feat_type]][[cell_a]][[cell_a]], ]

            predicted_values[[feat_type]][[feat_type]][[cell_a]] <- colSums(dimensions_cell_a) / k
        }
    }


    ######################## cross-modality prediction #########################

    if (verbose) message("Calculating cross-modality prediction")

    ## predicted feat_type1 modality2
    feat_type1 <- feat_types[1]
    feat_type2 <- feat_types[2]

    predicted_values[[feat_type1]][[feat_type2]] <- list()

    for (cell_a in cell_names) {
        dimensions_cell_a <- pca_list[[feat_type1]][kNN_list[[feat_type2]][[cell_a]][[cell_a]], ]

        predicted_values[[feat_type1]][[feat_type2]][[cell_a]] <- colSums(dimensions_cell_a) / k
    }

    ## predicted modality2 feat_type1
    predicted_values[[feat_type2]][[feat_type1]] <- list()

    for (cell_a in cell_names) {
        dimensions_cell_a <- pca_list[[feat_type2]][kNN_list[[feat_type1]][[cell_a]][[cell_a]], ]

        predicted_values[[feat_type2]][[feat_type1]][[cell_a]] <- colSums(dimensions_cell_a) / k
    }

    ###################### calculate jaccard similarities ######################

    if (verbose) message("Calculating Jaccard similarities")

    sNN_list <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        sNN_result <- createNearestNetwork(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            type = "sNN",
            dim_reduction_to_use = reduction_methods[i],
            dim_reduction_name = reduction_names[i],
            dimensions_to_use = 1:100,
            return_gobject = FALSE,
            minimum_shared = 1,
            k = 20
        )

        sNN_list[[feat_type]] <- igraph::as_data_frame(sNN_result)
    }

    if (verbose) message("Calculating kernel bandwidths")

    # cell-specific kernel bandwidth.
    sigma_i <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        sigma_i[[feat_type]] <- numeric()


        for (cell_a in cell_names) {
            ### 20 small jaccard values
            jaccard_values <- sNN_list[[feat_type]][sNN_list[[feat_type]]$from == cell_a, ]

            if (nrow(jaccard_values) == 20) {
                further_cell_cell_distances <- all_cell_distances[[feat_type]][[feat_type]][
                    cell_a, jaccard_values$to
                ]
            } else {
                further_cell_cell_distances <- tail(sort(all_cell_distances[[feat_type]][[feat_type]][
                    cell_a,
                ]), 20)
            }

            sigma_i[[feat_type]][cell_a] <- mean(further_cell_cell_distances)
            # cell-specific kernel bandwidth.
        }
    }

    ###################### cell-specific modality weights ######################

    if (verbose) message("Calculating modality weights")

    ## feat_type1 feat_type1
    theta_list <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        theta_list[[feat_type]] <- list()
        theta_list[[feat_type]][[feat_type]] <- list()

        for (cell_a in cell_names) {
            feat_type1_i <- pca_list[[feat_type]][cell_a, ] # profile of current cell

            d_feat_type1_i_modality2_predicted <- sqrt(sum((
                feat_type1_i - predicted_values[[feat_type]][[feat_type]][[cell_a]])^2))

            first_knn <- names(sort(cell_distances[[feat_type]][[feat_type]][[cell_a]]))[1]

            feat_type1_knn1 <- pca_list[[feat_type]][first_knn, ] # profile of the nearest neighbor

            d_feat_type1_i_feat_type1_knn1 <- sqrt(sum((
                feat_type1_i - feat_type1_knn1)^2))

            difference_distances <- d_feat_type1_i_modality2_predicted -
                d_feat_type1_i_feat_type1_knn1
            max_value <- max(c(difference_distances, 0))

            theta_list[[feat_type]][[feat_type]][[cell_a]] <- exp((
                -max_value) / (sigma_i[[feat_type]][cell_a] -
                d_feat_type1_i_feat_type1_knn1))
        }
    }


    ## feat_type1 modality2
    feat_type1 <- feat_types[1]
    feat_type2 <- feat_types[2]

    theta_list[[feat_type1]][[feat_type2]] <- list()

    for (cell_a in cell_names) {
        feat_type1_i <- pca_list[[feat_type1]][cell_a, ] # profile of current cell
        d_feat_type1_i_modality2_predicted <- sqrt(sum((
            feat_type1_i - predicted_values[[feat_type1]][[feat_type2]][[cell_a]])^2))
        first_knn <- names(sort(cell_distances[[feat_type1]][[feat_type1]][[cell_a]]))[1]

        feat_type1_knn1 <- pca_list[[feat_type1]][first_knn, ] # profile of the nearest neighbor
        d_feat_type1_i_feat_type1_knn1 <- sqrt(sum((
            feat_type1_i - feat_type1_knn1)^2))

        difference_distances <- d_feat_type1_i_modality2_predicted -
            d_feat_type1_i_feat_type1_knn1
        max_value <- max(c(difference_distances, 0))

        theta_list[[feat_type1]][[feat_type2]][[cell_a]] <- exp((
            -max_value) / (sigma_i[[feat_type1]][cell_a] -
            d_feat_type1_i_feat_type1_knn1))
    }


    ## modality2 feat_type1
    theta_list[[feat_type2]][[feat_type1]] <- list()

    for (cell_a in cell_names) {
        modality2_i <- pca_list[[feat_type2]][cell_a, ] # profile of current cell
        d_modality2_i_feat_type1_predicted <- sqrt(sum((
            modality2_i - predicted_values[[feat_type2]][[feat_type1]][[cell_a]])^2))
        first_knn <- names(sort(cell_distances[[feat_type2]][[feat_type2]][[cell_a]]))[1]

        modality2_knn1 <- pca_list[[feat_type2]][first_knn, ] # profile of the nearest neighbor
        d_modality2_i_modality2_knn1 <- sqrt(sum((
            modality2_i - modality2_knn1)^2))

        difference_distances <- d_modality2_i_feat_type1_predicted -
            d_modality2_i_modality2_knn1
        max_value <- max(c(difference_distances, 0))

        theta_list[[feat_type2]][[feat_type1]][[cell_a]] <- exp((
            -max_value) / (sigma_i[[feat_type2]][cell_a] -
            d_modality2_i_modality2_knn1))
    }


    ##################### ratio of affinities ##################################

    if (verbose) message("Calculating WNN")

    epsilon <- 10^-4

    ratio_list <- list()

    ## feat_type1
    ratio_list[[feat_type1]] <- list()

    for (cell_a in cell_names) {
        ratio_list[[feat_type1]][[cell_a]] <- theta_list[[feat_type1]][[feat_type1]][[cell_a]] /
            (theta_list[[feat_type1]][[feat_type2]][[cell_a]] + epsilon)
    }

    ## modality2
    ratio_list[[feat_type2]] <- list()

    for (cell_a in cell_names) {
        ratio_list[[feat_type2]][[cell_a]] <- theta_list[[feat_type2]][[feat_type2]][[cell_a]] /
            (theta_list[[feat_type2]][[feat_type1]][[cell_a]] + epsilon)
    }


    ########################### normalization ##################################

    if (verbose) message("Calculating WNN normalization")

    w_list <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        w_list[[feat_type]] <- rep(0, length(cell_names))
        names(w_list[[feat_type]]) <- cell_names

        for (cell_a in cell_names) {
            w_list[[feat_type]][cell_a] <- exp(ratio_list[[feat_type]][[cell_a]]) /
                (exp(ratio_list[[feat_type1]][[cell_a]]) + exp(ratio_list[[feat_type2]][[cell_a]]))
        }
    }

    ######################### Calculating a WNN graph ##########################

    if (verbose) message("Calculating WNN graph")

    theta_weighted <- matrix(rep(0, length(cell_names) * length(cell_names)),
        ncol = length(cell_names),
        nrow = length(cell_names)
    )

    colnames(theta_weighted) <- cell_names
    rownames(theta_weighted) <- cell_names

    kernelpower <- 1

    ## theta_feat_type1

    theta_cella_cellb <- list()

    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        theta_cella_cellb[[feat_type]] <- exp(-1 * (all_cell_distances[[feat_type]][[feat_type]] /
            sigma_i[[feat_type]])**kernelpower)
    }

    ## theta_weighted
    theta_weighted <- w_list[[feat_type1]] * theta_cella_cellb[[feat_type1]] +
        w_list[[feat_type2]] * theta_cella_cellb[[feat_type2]]


    # save theta_weighted

    if (verbose) message("Saving WNN results")


    ## set integrated feat_type and result name if not provided
    if (is.null(integrated_feat_type)) {
        integrated_feat_type <- paste0(feat_type1, "_", feat_type2)
    }

    if (is.null(matrix_result_name)) {
        matrix_result_name <- "theta_weighted_matrix"
    }

    gobject <- setMultiomics(
        gobject = gobject,
        x = theta_weighted,
        spat_unit = spat_unit,
        feat_type = integrated_feat_type,
        integration_method = "WNN",
        result_name = matrix_result_name,
        verbose = TRUE
    )

    # save feat_types weight
    for (i in seq(length(feat_types))) {
        feat_type <- feat_types[i]

        if (is.null(w_names[i])) {
            w_names[i] <- paste0("w_", feat_type)
        }

        gobject <- setMultiomics(
            gobject = gobject,
            x = w_list[feat_type],
            spat_unit = spat_unit,
            feat_type = integrated_feat_type,
            integration_method = "WNN",
            result_name = w_names[i],
            verbose = TRUE
        )
    }

    return(gobject)
}


# runIntegratedUMAP() of Giotto 4.3.6 (graph = "dense"), unchanged
.wnn_umap_dense <- function(
        gobject,
        spat_unit = "cell",
        feat_types = c("rna", "protein"),
        integrated_feat_type = NULL,
        integration_method = "WNN",
        matrix_result_name = "theta_weighted_matrix",
        k = 20,
        spread = 5,
        min_dist = 0.01,
        force = FALSE,
        seed = 1234,
        ...) {
    # validate feat_types
    for (feat_type in feat_types) {
        if (!feat_type %in% names(
            slot(gobject, "dimension_reduction")$cells[[spat_unit]]
        )) {
            stop(paste(feat_type, " and their dimension reductions must exist in the Giotto object"))
        }
    }

    feat_type1 <- feat_types[1]
    feat_type2 <- feat_types[2]

    if (is.null(integrated_feat_type)) {
        integrated_feat_type <- paste0(feat_type1, "_", feat_type2)
    }

    theta_weighted <- getMultiomics(gobject,
        spat_unit = spat_unit,
        feat_type = integrated_feat_type,
        integration_method = integration_method,
        result_name = matrix_result_name
    )

    theta_weighted[is.na(theta_weighted)] <- 0

    if (is.null(gobject@nn_network[[spat_unit]][[
        feat_type1
    ]]$kNN$integrated_kNN) || force == TRUE) {
        ################# Calculate integrated Nearest Neighbors ###############

        message("Calculating integrated Nearest Neighbors")

        cell_names <- colnames(theta_weighted)

        nn_network <- dbscan::kNN(x = theta_weighted, k = k, sort = TRUE)
        from <- to <- weight <- distance <- from_cell_ID <- to_cell_ID <-
            shared <- NULL
        nn_network_dt <- data.table::data.table(
            from = rep(
                seq_len(nrow(nn_network$id)),
                k
            ),
            to = as.vector(nn_network$id),
            weight = 1 / (1 + as.vector(nn_network$dist)),
            distance = as.vector(nn_network$dist)
        )
        nn_network_dt[, `:=`(from_cell_ID, cell_names[from])]
        nn_network_dt[, `:=`(to_cell_ID, cell_names[to])]
        all_index <- unique(
            x = c(nn_network_dt$from_cell_ID, nn_network_dt$to_cell_ID)
        )

        ################################ Create igraph #########################

        message("Creating igraph")

        nn_network_igraph <- igraph::graph_from_data_frame(
            nn_network_dt[, .(from_cell_ID, to_cell_ID, weight, distance)],
            directed = TRUE,
            vertices = all_index
        )

        ## store igraph
        nnNetObj <- createNearestNetObj(
            name = "integrated_kNN",
            nn_type = "kNN",
            network = nn_network_igraph,
            spat_unit = spat_unit,
            feat_type = feat_type1
        )

        gobject <- setGiotto(gobject, nnNetObj)

        ## store nn_network id
        gobject <- setMultiomics(
            gobject = gobject,
            x = nn_network$id,
            spat_unit = spat_unit,
            feat_type = integrated_feat_type,
            integration_method = "WNN",
            result_name = "integrated_kNN_id",
            verbose = TRUE
        )

        ## store nn_network dist
        gobject <- setMultiomics(
            gobject = gobject,
            x = nn_network$dist,
            spat_unit = spat_unit,
            feat_type = integrated_feat_type,
            integration_method = "WNN",
            result_name = "integrated_kNN_dist",
            verbose = TRUE
        )
    }

    ######################### Calculate integrated UMAP ########################

    message("Calculating integrated UMAP")

    nn_network_id <- getMultiomics(gobject,
        spat_unit = spat_unit,
        feat_type = integrated_feat_type,
        integration_method = integration_method,
        result_name = "integrated_kNN_id"
    )

    nn_network_dist <- getMultiomics(gobject,
        spat_unit = spat_unit,
        feat_type = integrated_feat_type,
        integration_method = integration_method,
        result_name = "integrated_kNN_dist"
    )


    #### using nn_network pre-calculation
    #
    # The stored graph comes from dbscan::kNN(), which REMOVES self-matches,
    # while uwot requires each observation to be its own first neighbor and
    # drops column 1 unconditionally when fitting the local connectivity
    # offset. Handing the stored matrices over unchanged therefore discarded
    # every cell's true nearest neighbor and fit rho against a shifted
    # distance set -- silently, because uwot validates neither the self
    # column nor the ordering. nnToUwot() adds the self entry back.
    #
    # This changes the integrated embedding. It was wrong before.
    GiottoUtils::local_seed(seed)
    integrated_umap <- uwot::umap(
        X = theta_weighted,
        n_neighbors = k,
        nn_method = GiottoClass::nnToUwot(
            list(id = nn_network_id, dist = nn_network_dist)
        ),
        spread = spread,
        min_dist = min_dist,
        ...
    )

    colnames(integrated_umap) <- c("Dim.1", "Dim.2")

    ## add umap
    for (feat_type in feat_types) {
        umap_object <- createDimObj(
            coordinates = integrated_umap,
            spat_unit = spat_unit,
            feat_type = feat_type,
            method = "umap",
            name = "integrated.umap"
        )

        gobject <- setDimReduction(
            gobject = gobject,
            x = umap_object,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction_method = "umap",
            name = "integrated.umap"
        )
    }

    return(gobject)
}
