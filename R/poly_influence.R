#' @title showPolygonSizeInfluence
#' @name showPolygonSizeInfluence
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param alt_spat_unit alternaitve spatial unit which represents resized
#' polygon data
#' @param feat_type feature type
#' @param clus_name name of cluster column in cell_metadata for given spat_unit
#' and alt_spat_unit, i.e. "kmeans"
#' @param return_plot logical. whether to return the plot object
#' @param verbose be verbose
#' @returns giotto object with altered cell_metadata
#' @details
#' Compares cell metadata from spat_unit-feat_type pairs as provided.
#'
#' New columns, resize_switch and cluster_interaction, will be created within
#' cell_metadata for spat_unit-feat_type.
#'
#' These new columns will describe if a given cell switched cluster number when
#' resized.
#' If the same amount of clusters exist for spat_unit-feat_type and
#' alt_spat_unit-feat_type, then clusters are determined to be
#' corresponding based on % overlap in cell_IDs in each cluster.
#'
#' Otherwise, multiple clusters from the spatial unit feature type pair are
#' condensed to align with the smaller number of clusters and ensure overlap.
#'
#' @export
showPolygonSizeInfluence <- function(
        gobject = NULL,
        spat_unit = NULL,
        alt_spat_unit = NULL,
        feat_type = NULL,
        clus_name = "kmeans",
        return_plot = FALSE,
        verbose = FALSE) {
    # NSE vars
    cell_ID <- total_expr <- cluster_interactions <- N <- resize_switch <- NULL

    # Guards
    if (!c("giotto") %in% class(gobject)) {
        stop(wrap_txt("Please provide a valid Giotto Object.", errWidth = TRUE))
    }

    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    if (!alt_spat_unit %in% names(gobject@expression)) {
        stop(wrap_txt(paste0(
            "Alternative spatial unit ", alt_spat_unit,
            " not found. Please ensure it exists."
        ), errWidth = TRUE))
    }

    meta_cols <- names(getCellMetadata(
        gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    ))

    if (!clus_name %in% meta_cols) {
        stop(wrap_txt(
            paste0(
                "Cluster name ", clus_name,
                " not found within cell metadata. Please ensure it exists."
            ),
            errWidth = TRUE
        ))
    }


    if (c("cluster_interactions") %in% meta_cols) {
        warning((wrap_txt(paste0(
            "Switch interactions already found within
        cell_metadata for spat_unit feat_type pair:`", spat_unit, "-",
            feat_type, "`. They will be overwritten."
        ), errWidth = TRUE)))
    }
    ## Compare clustering results between cell and smallcell data #######
    # ----------------------------------------------------------------- #
    cell_meta <- pDataDT(gobject, spat_unit = spat_unit)
    cell_meta_new <- pDataDT(gobject, spat_unit = alt_spat_unit)

    new_clus_table <- cell_meta_new[, .(cell_ID, total_expr)]
    new_clus_table[[clus_name]] <- cell_meta_new[[clus_name]]

    cell_meta <- merge.data.table(cell_meta, new_clus_table, by = "cell_ID")

    cell_meta[, cluster_interactions := paste0(cell_meta[[
        paste0(clus_name, ".x")
    ]], "-", cell_meta[[paste0(clus_name, ".y")]])]
    switches2 <- cell_meta[, .N, by = "cluster_interactions"]
    setorder(switches2, N)

    num_orig <- sort(unique(cell_meta[[paste0(clus_name, ".x")]]))
    num_new <- sort(unique(cell_meta[[paste0(clus_name, ".y")]]))

    equal_len <- TRUE
    if (length(num_orig) != length(num_new)) equal_len <- FALSE

    switch_strs <- c() # scope
    if (!equal_len) {
        #####
        switch_strs <- .determine_switch_string_unequal(
            num_orig = num_orig,
            num_new = num_new
        )
        #####
    } else {
        cmeta <- pDataDT(gobject, spat_unit = spat_unit)
        cmeta_new <- pDataDT(gobject, spat_unit = alt_spat_unit)

        switch_strs <- .determine_switch_string_equal(
            cell_meta = cmeta,
            cell_meta_new = cmeta_new,
            clus_name = clus_name
        )
    }

    cell_meta[, resize_switch := ifelse(
        cluster_interactions %in% switch_strs, "same", "switch"
    )]
    gobject <- addCellMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        new_metadata = cell_meta[
            , .(cell_ID, resize_switch, cluster_interactions)
        ],
        by_column = TRUE,
        column_cell_ID = "cell_ID"
    )

    poly_plot <- spatInSituPlotPoints(
        gobject = gobject,
        spat_unit = spat_unit,
        polygon_feat_type = spat_unit,
        show_polygon = TRUE,
        feat_type = feat_type,
        feats = NULL,
        polygon_fill = "resize_switch",
        polygon_fill_as_factor = TRUE,
        polygon_line_size = 0.1,
        polygon_color = "white",
        coord_fix_ratio = 1,
        polygon_fill_code = c(switch = "red", same = "gray"),
        return_plot = return_plot
    )

    num_cells_switched <- sum(
        getCellMetadata(gobject)$resize_switch == "switch"
    )
    num_cells_same <- sum(getCellMetadata(gobject)$resize_switch == "same")
    if (verbose) print(paste0(num_cells_switched, " cells switched clusters."))
    if (verbose) {
        print(paste0(num_cells_same, " cells remained in the same cluster."))
    }

    if (return_plot) {
        return(poly_plot)
    }

    return(gobject)
}

#' @title Determine switch string equal
#' @name .determine_switch_string_equal
#' @param cell_meta cell_metadata from the original spatial_unit
#' @param cell_meta_new cell_metadata from the resized spatial unit
#' @param clus_name name of the cluster type, likely "kmeans"
#' @returns switch_str, a vector of corresponding cluster numbers in strings
#' @details creates a string in the format c("x_1-y_1", "x_2-y_2"..."x_n, y_n")
#' Where:
#'    x_n is a cluster number from the original spatial unit
#'    y_m is a cluster number from the resized spatial unit
#'    n is the number of clusters
#'
#' Clusters are determined to be corresponding based on % overlap in cell_IDs
#' in each cluster.
#'
#' @keywords internal
.determine_switch_string_equal <- function(
        cell_meta = NULL,
        cell_meta_new = NULL,
        clus_name = NULL) {
    k_clusters <- sort(unique(cell_meta[[clus_name]]))
    num_clusters <- k_clusters[length(k_clusters)]

    k_match_clusters <- seq_len(num_clusters)
    switch_strs <- c()
    for (i in seq_len(num_clusters)) {
        thresh <- 0
        clus_match <- NULL
        for (j in seq_len(num_clusters)) {
            c_df <- cell_meta[cell_meta[[clus_name]] == i]$cell_ID
            nc_df <- cell_meta_new[cell_meta_new[[clus_name]] == j]$cell_ID

            overlap <- sum(c_df %in% nc_df / length(c_df))
            if (overlap > thresh) {
                thresh <- overlap
                clus_match <- j
            }
        }

        k_match_clusters[i] <- clus_match
    }

    for (idx in seq_len(num_clusters)) {
        p1 <- k_clusters[[idx]]
        p2 <- k_match_clusters[[idx]]
        switch_strs <- c(switch_strs, paste0(p1, "-", p2))
    }

    return(switch_strs)
}

#' @title Determine switch string unequal
#' @name .determine_switch_string_unequal
#' @param num_orig sorted vector of cluster numbers in the original metadata
#' @param num_new sorted vector of cluster numbers in the new, resized metadata
#' @returns switch_str, a vector of corresponding cluster numbers in strings
#' @details determines how to create a string in the format
#' c("x_1-y_1", "x_2-y_2"..."x_n, y_m")
#' Where:
#'    x_n is a cluster number from the original spatial unit
#'    y_m is a cluster number from the resized spatial unit
#'    n is the number of clusters in the original spatial unit
#'    m is the number of clusters in the new spatial unit
#'
#' Essentially determines iteration order for .create_switch_string_unequal()
#'
#' @keywords internal
.determine_switch_string_unequal <- function(
        num_orig = NULL,
        num_new = NULL) {
    switch_strs <- c()

    orig_first <- TRUE

    if (length(num_orig) < length(num_new)) orig_first <- FALSE

    if (orig_first) {
        switch_strs <- .create_switch_string_unequal(
            num_first = num_orig,
            num_second = num_new,
            switch_strs = switch_strs
        )
        return(switch_strs)
    }

    switch_strs <- .create_switch_string_unequal(
        num_first = num_new,
        num_second = num_orig,
        switch_strs = switch_strs
    )

    return(switch_strs)
}

#' @title Create switch string unequal
#' @name .create_switch_string_unequal
#' @param num_first sorted vector of cluster numbers in the outer for loop
#' @param num_second sorted vector of cluster numbers in the inner for loop
#' @returns switch_str, a vector of corresponding cluster numbers in strings
#' @details creates a string in the format c("x_1-y_1", "x_2-y_2"..."x_n, y_m")
#' Where:
#'    x_n is a cluster number from the original spatial unit
#'    y_m is a cluster number from the resized spatial unit
#'    n is the number of clusters in the original spatial unit
#'    m is the number of clusters in the new spatial unit
#' @keywords internal
.create_switch_string_unequal <- function(
        num_first = NULL,
        num_second = NULL,
        switch_strs = NULL) {
    for (o in num_first) {
        for (n in num_second) {
            if (as.integer(o) == as.integer(n)) {
                switch_strs <- c(
                    switch_strs, paste0(as.character(o), "-", as.character(n))
                )
            }
            if (o > n && n == num_second[length(num_second)]) {
                switch_strs <- c(
                    switch_strs, paste0(as.character(o), "-", as.character(n))
                )
            }
        }
    }

    switch_strs <- unique(switch_strs)

    return(switch_strs)
}

#' @title showCellProportionSwitchedPie
#' @name showCellProportionSwitchedPie
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @returns ggplot
#' @details Creates a pie chart showing how many cells switched clusters after
#' annotation resizing.
#' The function showPolygonSizeInfluence() must have been run on the Giotto
#' Object for this function to run.
#' @export
showCellProportionSwitchedPie <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL) {
    # NSE vars
    cluster_status <- num_cells <- resize_switch <- perc <- ypos <- NULL

    # Guards
    if (!"giotto" %in% class(gobject)) {
        stop(wrap_txt("Please provide a valid Giotto Object.", errWidth = TRUE))
    }

    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )


    # Extract cell metadata
    cmeta <- getCellMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )

    if (!c("resize_switch") %in% names(cmeta)) {
        stop(wrap_txt("Column 'resize_switch' not found in cell metadata.
            Ensure showPolygonSizeInfluence() has been run.", errWidth = TRUE))
    }

    plotdf <- data.table::data.table()
    plotdf[, cluster_status := c("switch", "same")]
    plotdf[, num_cells := c(sum(cmeta[
        , resize_switch == "switch"
    ]), sum(cmeta[, resize_switch == "same"]))]

    per_switch <- plotdf$num_cells[[1]] / sum(plotdf$num_cells) * 100
    per_same <- plotdf$num_cells[[2]] / sum(plotdf$num_cells) * 100

    y_switch <- cumsum(per_switch) - 0.5 * per_switch
    y_same <- cumsum(per_same + per_switch) - 0.5 * per_same


    plotdf[, perc := c(per_switch, per_same)]
    plotdf[, ypos := c(y_switch, y_same)]

    print(plotdf)

    ggplot(
        as.data.frame(plotdf), aes(x = "", y = perc, fill = cluster_status)
    ) +
        coord_polar("y", start = 0) +
        geom_bar(stat = "identity", width = 1) +
        theme_void() +
        geom_text(aes(y = ypos, label = num_cells), color = "white", size = 6)
}

#' @title showCellProportionSwitchedSanKey
#' @name showCellProportionSwitchedSanKey
#' @param gobject giotto object which contains metadata for both spat_unit and
#' alt_spat_unit
#' @param spat_unit spatial unit
#' @param alt_spat_unit alternative spatial unit which stores data after
#' resizing annotations
#' @param feat_type feature type
#' @returns D3 JavaScript Sankey diagram
#' @details Creates a Sankey Diagram to illustrate cluster switching behavior.
#' Currently only supports displaying cluster switching for kmeans clusters.
#' @export
showCellProportionSwitchedSanKey <- function(
        gobject = NULL,
        spat_unit = NULL,
        alt_spat_unit = NULL,
        feat_type = NULL) {
    # NSE vars
    kmeans_small <- cell_ID <- NULL

    # Guards
    if (!"giotto" %in% class(gobject)) {
        stop(wrap_txt("Please provide a valid Giotto Object.", errWidth = TRUE))
    }

    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )
    if (!alt_spat_unit %in% names(gobject@expression)) {
        stop(wrap_txt(paste0(
            "Alternative spatial unit ", alt_spat_unit,
            " not found. Please ensure it exists."
        ), errWidth = TRUE))
    }

    package_check("networkD3")

    # Extract cell metadata
    cmeta <- getCellMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )

    if (!c("resize_switch") %in% names(cmeta)) {
        stop(wrap_txt("Column 'resize_switch' not found in cell metadata.
            Ensure showPolygonSizeInfluence() has been run.", errWidth = TRUE))
    }

    if (!c("kmeans") %in% names(cmeta)) {
        stop("This function has only been implemented for k-means clusters.")
    }


    small_cmeta <- getCellMetadata(
        gobject = gobject,
        spat_unit = alt_spat_unit,
        feat_type = feat_type,
        output = "data.table"
    )

    if (!c("kmeans") %in% names(small_cmeta)) {
        stop("This function has only been implemented for k-means clusters.")
    }

    small_cmeta_clus <- small_cmeta[, .(cell_ID, kmeans)]
    small_cmeta_clus$kmeans_small <- small_cmeta_clus$kmeans
    small_cmeta_clus$kmeans <- NULL

    merged_cmeta <- data.table::merge.data.table(
        cmeta, small_cmeta_clus,
        by.x = "cell_ID", by.y = "cell_ID"
    )


    k1 <- unique(merged_cmeta$kmeans)
    k2 <- unique(merged_cmeta$kmeans_small)

    fdt <- data.table::data.table()
    c_k1 <- c()
    c_k2 <- c()

    flen <- length(k1) * length(k2)
    idx1 <- 1
    idx2 <- 1

    for (i in seq_len(flen)) {
        c_k1[i] <- k1[idx1] - 1 # java zero-index
        c_k2[i] <- k2[idx2] - 1 # java zero-index

        if (i %% length(k1) == 0) idx1 <- idx1 + 1
        idx2 <- idx2 + 1
        if (idx2 > length(k2)) idx2 <- 1
    }

    num_occ <- c()

    for (i in seq_len(flen)) {
        num_occ[i] <- dim(na.omit(merged_cmeta[kmeans == (c_k1[i] + 1)][
            merged_cmeta[kmeans_small == (c_k2[i] + 1)]
        ]))[[1]]
    }

    fdt[, "k1"] <- c_k1
    fdt[, "k2"] <- c_k2 + 7
    fdt[, "value"] <- num_occ
    fdt

    label_dt <- data.table::data.table()
    label_dt[, "name"] <- c(paste0("original_", as.character(sort(k1))), paste0(
        "resized_", as.character(sort(k2))
    ))
    label_dt

    master <- list(fdt, label_dt)
    names(master) <- c("links", "nodes")

    networkD3::sankeyNetwork(
        Links = master$links,
        Nodes = master$nodes,
        Source = "k1",
        Target = "k2",
        Value = "value",
        NodeID = "name",
        units = "TWh",
        fontSize = 12,
        nodeWidth = 30
    )
}
