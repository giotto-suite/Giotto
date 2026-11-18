#' Spatial domains with BANKSY
#'
#' @description
#' Runs BANKSY (Singhal et al. 2024, *Nat Genet* 56:431) through the
#' Bioconductor package \pkg{Banksy} on a stored dimension reduction: each
#' cell's coordinates are joined by the mean of its spatial neighbors'
#' (and, with `compute_agf = TRUE`, their azimuthal gradient), weighted by
#' `lambda`, and reduced by PCA; the BANKSY PCs are then clustered with
#' Leiden into domains. The defaults are the BANKSY authors' settings for
#' finding domains. With a joint embedding such as `runBlockPCA()`'s, the
#' domains use RNA and protein together.
#'
#' Only the embedding (cells x dimensions) and the spatial locations are
#' passed to Banksy, never the expression, so GiottoDisk-backed objects work
#' as in-memory ones.
#'
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type of the dimension reduction; the domains are
#'   added to its cell metadata
#' @param dim_reduction_to_use reduction method of the embedding, e.g.
#'   `"pca"` or `"blockpca"`
#' @param dim_reduction_name name of the embedding. `NULL` takes the default
#'   name of `dim_reduction_to_use` for `feat_type`.
#' @param dimensions_to_use dimensions (indices) of the embedding to use.
#'   `NULL` (default) uses all; indices beyond those stored are dropped.
#' @param spat_loc_name name of the spatial locations
#' @param lambda weight of the neighborhood in `[0, 1]`: about 0.8 finds
#'   spatial domains, about 0.2 refines cell types with their surroundings
#' @param k_geom number of spatial neighbors for the neighborhood mean and,
#'   with `compute_agf = TRUE`, for the azimuthal gradient
#' @param compute_agf logical. Also use the azimuthal Gabor filter (gradient)
#'   of the neighborhood
#' @param npcs number of PCs of the BANKSY matrix used for clustering
#' @param cluster_method `"giotto"` (default): a shared nearest neighbor
#'   network of the BANKSY PCs ([createNearestNetwork()], stored as `name`)
#'   and [doLeidenCluster()] run to convergence. `"banksy"`:
#'   `Banksy::clusterBanksy()`, BANKSY's own graph and Leiden; on 146,000 cells
#'   it takes about 5 minutes on one thread and about 10 GB for its graph.
#' @param k_neighbors neighbors of the clustering network
#' @param nn_engine `"giotto"` clustering: neighbor search of
#'   [createNearestNetwork()], exact (`"dbscan"`) or approximate (`"hnsw"`,
#'   faster on whole slides)
#' @param resolution Leiden resolution
#' @param name name of the cell metadata column with the domains, and of the
#'   stored BANKSY PCs (reduction method `"banksy"`)
#' @param seed seed
#' @param verbose be verbose
#' @returns The Giotto object with the domains in the cell metadata of
#'   `feat_type` (column `name`) and the BANKSY PCs as a dimension reduction
#'   (method `"banksy"`, name `name`); with `cluster_method = "giotto"` also
#'   the sNN network `name`.
#' @seealso `runBlockPCA()` for a joint RNA + protein embedding.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' g <- runBanksy(g, dimensions_to_use = 1:20, npcs = 10, k_geom = 6,
#'     resolution = 0.5)
#' table(pDataDT(g)$banksy)
#' @md
#' @export
runBanksy <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        dim_reduction_to_use = "pca",
        dim_reduction_name = NULL,
        dimensions_to_use = NULL,
        spat_loc_name = "raw",
        lambda = 0.8,
        k_geom = 30,
        compute_agf = FALSE,
        npcs = 20,
        cluster_method = c("giotto", "banksy"),
        k_neighbors = 50,
        nn_engine = c("dbscan", "hnsw"),
        resolution = 0.8,
        name = "banksy",
        seed = 1234,
        verbose = NULL) {
    checkmate::assert_class(gobject, "giotto")
    checkmate::assert_string(dim_reduction_to_use)
    checkmate::assert_string(dim_reduction_name, null.ok = TRUE)
    checkmate::assert_integerish(dimensions_to_use, lower = 1,
        any.missing = FALSE, null.ok = TRUE)
    checkmate::assert_string(spat_loc_name)
    checkmate::assert_number(lambda, lower = 0, upper = 1)
    checkmate::assert_integerish(k_geom, lower = 1, any.missing = FALSE,
        min.len = 1L, max.len = 2L)
    checkmate::assert_flag(compute_agf)
    checkmate::assert_count(npcs, positive = TRUE)
    cluster_method <- match.arg(cluster_method)
    checkmate::assert_count(k_neighbors, positive = TRUE)
    nn_engine <- match.arg(nn_engine)
    checkmate::assert_number(resolution, lower = 0)
    checkmate::assert_string(name)
    package_check(c("SpatialExperiment", "Banksy"),
        repository = c("Bioc", "Bioc"))
    GiottoUtils::local_seed(seed)

    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    emb <- getDimReduction(gobject,
        spat_unit = spat_unit, feat_type = feat_type, reduction = "cells",
        reduction_method = dim_reduction_to_use, name = dim_reduction_name,
        output = "matrix"
    )
    if (!is.null(dimensions_to_use)) {
        dims <- dimensions_to_use[dimensions_to_use <= ncol(emb)]
        if (!length(dims)) {
            stop("[runBanksy] no stored dimension is selected",
                call. = FALSE)
        }
        emb <- emb[, dims, drop = FALSE]
    }
    if (any(!is.finite(emb))) {
        stop("[runBanksy] the embedding has NA or infinite values",
            call. = FALSE)
    }
    npcs <- min(npcs, ncol(emb) * (1L + length(k_geom) * compute_agf))
    sl <- getSpatialLocations(gobject,
        spat_unit = spat_unit, name = spat_loc_name, output = "data.table"
    )
    sl <- sl[match(rownames(emb), sl$cell_ID)]
    if (anyNA(sl$cell_ID)) {
        stop("[runBanksy] some cells of the embedding have no spatial ",
            "location", call. = FALSE)
    }

    vmsg(.v = verbose, "runBanksy: neighborhood matrices")
    se <- SpatialExperiment::SpatialExperiment(
        assays = list(embedding = t(emb)),
        spatialCoords = as.matrix(sl[, c("sdimx", "sdimy")])
    )
    se <- Banksy::computeBanksy(se, assay_name = "embedding",
        compute_agf = compute_agf, k_geom = k_geom, verbose = FALSE)
    vmsg(.v = verbose, "runBanksy: PCA")
    se <- Banksy::runBanksyPCA(se, use_agf = compute_agf, lambda = lambda,
        npcs = npcs, seed = seed)
    pcs <- SingleCellExperiment::reducedDim(se,
        grep("^PCA_", SingleCellExperiment::reducedDimNames(se),
            value = TRUE)[1L])
    dimnames(pcs) <- list(rownames(emb), paste0("Dim.", seq_len(ncol(pcs))))

    dim_obj <- createDimObj(
        coordinates = pcs, name = name, spat_unit = spat_unit,
        feat_type = feat_type, method = "banksy", reduction = "cells",
        misc = list(input = list(dim_reduction_to_use = dim_reduction_to_use,
            dim_reduction_name = dim_reduction_name,
            dimensions_to_use = colnames(emb)), lambda = lambda,
            k_geom = k_geom, compute_agf = compute_agf,
            cluster_method = cluster_method, k_neighbors = k_neighbors,
            resolution = resolution)
    )
    gobject <- setDimReduction(gobject, x = dim_obj, spat_unit = spat_unit,
        feat_type = feat_type, reduction = "cells",
        reduction_method = "banksy", name = name, verbose = FALSE)

    vmsg(.v = verbose, "runBanksy: Leiden")
    if (cluster_method == "giotto") {
        gobject <- createNearestNetwork(gobject, spat_unit = spat_unit,
            feat_type = feat_type, type = "sNN",
            dim_reduction_to_use = "banksy", dim_reduction_name = name,
            dimensions_to_use = seq_len(ncol(pcs)), k = k_neighbors,
            engine = nn_engine, keep_knn = FALSE, name = name,
            verbose = FALSE)
        return(doLeidenCluster(gobject, spat_unit = spat_unit,
            feat_type = feat_type, name = name, nn_network_to_use = "sNN",
            network_name = name, resolution = resolution, n_iterations = -1,
            set_seed = TRUE, seed_number = seed))
    }
    se <- Banksy::clusterBanksy(se, use_agf = compute_agf, lambda = lambda,
        npcs = npcs, k_neighbors = k_neighbors, resolution = resolution,
        seed = seed)
    cd <- SummarizedExperiment::colData(se)
    new_meta <- data.table::data.table(cell_ID = rownames(emb))
    new_meta[[name]] <- as.character(cd[[grep("^clust_", colnames(cd))[1L]]])
    addCellMetadata(gobject, spat_unit = spat_unit, feat_type = feat_type,
        new_metadata = new_meta, by_column = TRUE, column_cell_ID = "cell_ID")
}
