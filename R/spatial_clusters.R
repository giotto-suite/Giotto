#' @title Split cluster annotations based on a spatial network
#' @name spatialSplitCluster
#' @inheritParams data_access_params
#' @param spatial_network_name character. Name of spatial network to use
#' @param cluster_col character. Column in metadata containing original
#' clustering
#' @param split_clus_name character. Name to assign the split cluster results
#' @param include_all_ids logical. Include all ids, including vertex ids not
#' found in the spatial network
#' @param return_gobject logical. Return giotto object
#' @returns giotto object with cluster annotations
#' @examples
#' g <- GiottoData::loadGiottoMini("vizgen")
#' activeSpatUnit(g) <- "aggregate"
#' spatPlot2D(g, cell_color = "leiden_clus")
#'
#' g <- spatialSplitCluster(g,
#'     cluster_col = "leiden_clus",
#'     split_clus_name = "new"
#' )
#' # don't show legend since there are too many categories generated
#' spatPlot2D(g, cell_color = "new", show_legend = FALSE)
#' @export
spatialSplitCluster <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        spatial_network_name = "Delaunay_network",
        cluster_col,
        split_clus_name = paste0(cluster_col, "_split"),
        missing_id_name = "not_connected",
        return_gobject = TRUE) {
    # NSE vars
    cell_ID <- NULL

    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    cell_meta <- getCellMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = FALSE
    )

    sn <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        output = "spatialNetworkObj",
        copy_obj = FALSE,
        verbose = FALSE,
    )

    clus_info <- cell_meta[, c("cell_ID", cluster_col), with = FALSE]
    # subset to needed cols
    g <- GiottoClass::spat_net_to_igraph(sn)
    # convert spatialNetworkObject to igraph

    # assign cluster info to igraph nodes
    clus_values <- clus_info[
        match(igraph::V(g)$name, cell_ID), get(cluster_col)
    ]
    igraph::V(g)$cluster <- clus_values

    # split cluster by spatial igraph
    g <- .igraph_remove_hetero_edges(
        g = g,
        clus_attr = "cluster"
    )

    # get new clusterings
    # spatially unconnected nodes (if any) will always be returned as 0
    all_ids <- unique(cell_meta$cell_ID)
    new_clus_dt <- .igraph_vertex_membership(
        g = g,
        clus_name = split_clus_name,
        all_ids <- all_ids
    )

    if (isTRUE(return_gobject)) {
        gobject <- addCellMetadata(
            gobject,
            spat_unit = spat_unit,
            new_metadata = new_clus_dt,
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )
        return(gobject)
    } else {
        new_clus_dt
    }
}





#' @title Split cluster annotations based on a spatial network
#' @name identifyTMAcores
#' @inheritParams data_access_params
#' @param spatial_network_name character. Name of spatial network to use
#' @param core_id_name metadata column name for the core information
#' @param id_fmt character. [sprintf] formatting to use for core ids
#' @param include_all_ids logical. Include all ids, including vertex ids not
#' found in the spatial network
#' @param missing_id_name character. Name for nodes that are not connected to
#' a core.
#' @param min_nodes numeric. Minimal number of nodes to not be considered
#' an unconnected group.
#' @param join_split_cores logical. Attempt to repair core IDs when a core
#' is split down the middle and detected as two different cores.
#' @param join_tolerance numeric. Max ratio allowed relative to previous max
#' core convex hull area when determining if a pair of cores should be joined.
#' @param return_gobject logical. Return giotto object
#' @returns cluster annotations
#' @export
identifyTMAcores <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        spatial_network_name = "Delaunay_network",
        core_id_name = "core_id",
        id_fmt = "%d",
        include_all_ids = TRUE,
        missing_id_name = "not_connected",
        min_nodes = 5,
        join_split_cores = TRUE,
        join_tolerance = 1.2,
        return_gobject = TRUE) {
    # NSE vars
    cell_ID <- NULL

    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # get data
    cell_meta <- getCellMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = FALSE
    )

    sn <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        output = "spatialNetworkObj",
        copy_obj = FALSE,
        verbose = FALSE,
    )

    g <- GiottoClass::spat_net_to_igraph(sn)
    # convert spatialNetworkObject to igraph

    # get new clusterings as initial indices
    # these indices may need repairs and updates to be finalized
    ivm_params <- list(
        g = g, clus_name = "init_idx"
    )
    if (isTRUE(include_all_ids)) {
        # include all cell IDs
        all_ids <- unique(cell_meta$cell_ID)
        ivm_params$all_ids <- all_ids
    } else {
        # only IDs present in graph
        ivm_params$all_ids <- NULL
    }

    new_clus_dt <- do.call(.igraph_vertex_membership, args = ivm_params)
    # connected nodes
    con <- new_clus_dt[init_idx > 0]
    # spatially disconnected observations (not connected to a group of nodes)
    dcon <- new_clus_dt[init_idx == 0]

    # min nodes filter
    con_nodes <- con[, .N, by = init_idx]
    small_con_idx <- con_nodes[N < min_nodes, init_idx]
    # shift filtered values to dcon (disconnected)
    con[init_idx %in% small_con_idx, init_idx := 0]
    dcon <- rbind(dcon, con[init_idx == 0])
    con <- con[init_idx != 0]

    # fix split cores
    if (join_split_cores) {
        sl <- getSpatialLocations(gobject, spat_unit = spat_unit)
        con_init_idx_uniq <- sort(unique(con$init_idx))

        areas <- vapply(
            FUN.VALUE = numeric(1L), con_init_idx_uniq, function(core_id) {
                sl[con[init_idx == core_id, cell_ID]] |>
                    convHull() |>
                    area()
            }
        )
        max_area <- max(areas)

        # find ext of cores
        # iterate through angles to catch cases where extents do not
        # bridge across split.
        ovlp_reps <- lapply(c(0, 22.5, 45), function(rangle) {
            sl_rot <- spin(sl, rangle)

            # get ext poly of rotated cores
            epoly_list <- lapply(con_init_idx_uniq, function(core_id) {
                sl_rot[con[init_idx == core_id, cell_ID]] |>
                    ext() |>
                    as.polygons()
            })
            poly <- do.call(rbind, epoly_list)

            # test for overlaps
            ovlp <- relate(poly, relation = "overlaps", pairs = TRUE) |>
                # determine sorted pairs of overlaps
                apply(MARGIN = 2, sort) |>
                t()
            return(ovlp)
        })
        # combine test reps
        ovlps <- do.call(rbind, ovlp_reps) |>
            unique()

        # update ids based on test
        for (pair_i in nrow(ovlps)) {
            idx_1 <- ovlps[pair_i, 1L]
            idx_2 <- ovlps[pair_i, 2L]
            # ignore hits from two full cores
            # combined area of IDs to join cannot be greater than join_tolerance of max_area
            if ((areas[[idx_1]] + areas[[idx_2]]) >
                (join_tolerance * max_area)) {
                next
            }

            con[init_idx == idx_2, init_idx := idx_1]
        }
    }

    # apply core_id_name
    con[, (core_id_name) := sprintf(id_fmt, init_idx)]
    dcon[, (core_id_name) := missing_id_name]

    res <- rbind(
        con[, .SD, .SDcols = c("cell_ID", core_id_name)],
        dcon[, .SD, .SDcols = c("cell_ID", core_id_name)]
    )

    if (isTRUE(return_gobject)) {
        gobject <- addCellMetadata(
            gobject,
            spat_unit = spat_unit,
            new_metadata = res,
            by_column = TRUE,
            column_cell_ID = "cell_ID"
        )
        return(gobject)
    } else {
        new_clus_dt
    }
}












# internals ####

#' @title Remove hetero edges from igraph
#' @name .igraph_remove_hetero_edges
#' @description
#' Given an igraph `g` and set of node attributes `clus_att` that encode
#' different spatial clusters, remove edges that connect non-similar nodes.
#' This can be used when data is already clustered, but these clusters should
#' be further broken up based on whether they are spatially touching.
#' @param g igraph object
#' @param clus_attr character. A categorical node attribute
#' @returns igraph
#' @noRd
#' @keywords internal
.igraph_remove_hetero_edges <- function(g, clus_attr) {
    clus_attr_values <- igraph::vertex_attr(g, name = clus_attr)

    for (n in unique(clus_attr_values)) {
        # find all vertices of the attribute
        nv <- igraph::V(g)$name[clus_attr_values == n]

        # find edges that include these vertices
        n_all_edges <- igraph::E(g)[.inc(igraph::V(g)[nv])] %>%
            igraph::as_ids()

        # find edges associated with only these vertices
        n_internal_edges <- igraph::E(g)[nv %--% nv] %>%
            igraph::as_ids()

        het_edges <- n_all_edges[!n_all_edges %in% n_internal_edges]

        g <- igraph::delete_edges(g, edges = het_edges)
    }

    g
}




#' @title igraph vertex membership
#' @name .igraph_vertex_membership
#' @description
#' Get which weakly connected set of vertices each vertex is part of
#' @param g igraph
#' @param clus_name character. name to assign column of clustering info
#' @param all_ids (optional) character vector with all ids
#' @returns `data.table` with two columns. 1st is "cell_ID", second is named with
#' `clus_name` and is of type `numeric`
#' @keywords internal
#' @noRd
.igraph_vertex_membership <- function(g,
    clus_name,
    all_ids = NULL) {
    # get membership
    membership <- igraph::components(g)$membership %>%
        data.table::as.data.table(keep.rownames = TRUE)
    data.table::setnames(membership, c("cell_ID", clus_name))

    # add vertices that were missing from g back
    if (!is.null(all_ids)) {
        missing_ids <- all_ids[!all_ids %in% igraph::V(g)$name]
        missing_membership <- data.table::data.table(
            "cell_ID" = missing_ids,
            "cluster_name" = 0
        )
        data.table::setnames(missing_membership, c("cell_ID", clus_name))
        membership <- data.table::rbindlist(
            list(membership, missing_membership)
        )
    }

    return(membership)
}
