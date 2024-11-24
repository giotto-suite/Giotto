#' @title cellProximityBarplot
#' @name cellProximityBarplot
#' @description Create barplot from cell-cell proximity scores
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param CPscore CPscore, output from cellProximityEnrichment()
#' @param min_orig_ints filter on minimum original cell-cell interactions
#' @param min_sim_ints filter on minimum simulated cell-cell interactions
#' @param p_val p-value
#' @returns ggplot barplot
#' @details This function creates a barplot that shows the  spatial proximity
#'  enrichment or depletion of cell type pairs.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' cellProximityBarplot(
#'     gobject = g,
#'     CPscore = cellProximityEnrichment(g, cluster_column = "leiden_clus")
#' )
#' @export
cellProximityBarplot <- function(
        gobject,
        CPscore,
        min_orig_ints = 5,
        min_sim_ints = 5,
        p_val = 0.05,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "cellProximityBarplot") {
    table_mean_results_dc <- CPscore$enrichm_res

    ## filter to remove low number of cell-cell proximity interactions ##

    # data.table variables
    original <- simulations <- p_higher_orig <- p_lower_orig <- enrichm <-
        type_int <- unified_int <- NULL

    table_mean_results_dc_filter <- table_mean_results_dc[
        original >= min_orig_ints & simulations >= min_sim_ints,
    ]
    table_mean_results_dc_filter <- table_mean_results_dc_filter[
        p_higher_orig <= p_val | p_lower_orig <= p_val,
    ]

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_bar(
        data = table_mean_results_dc_filter,
        ggplot2::aes(x = unified_int, y = enrichm, fill = type_int),
        stat = "identity", show.legend = FALSE
    )
    pl <- pl + ggplot2::coord_flip()
    pl <- pl + ggplot2::theme_bw()
    pl <- pl + ggplot2::labs(y = "enrichment/depletion")
    pl

    bpl <- ggplot2::ggplot()
    bpl <- bpl + ggplot2::geom_bar(
        data = table_mean_results_dc_filter,
        ggplot2::aes(x = unified_int, y = original, fill = type_int),
        stat = "identity", show.legend = TRUE
    )
    bpl <- bpl + ggplot2::coord_flip()
    bpl <- bpl + ggplot2::theme_bw() + ggplot2::theme(
        axis.text.y = element_blank()
    )
    bpl <- bpl + ggplot2::labs(y = "# of interactions")
    bpl

    combo_plot <- plot_grid(
        pl, bpl,
        ncol = 2, rel_heights = c(1),
        rel_widths = c(3, 1.5), align = "h"
    )

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = combo_plot,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}

#' @title cellProximityHeatmap
#' @name cellProximityHeatmap
#' @description Create heatmap from cell-cell proximity scores
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param CPscore CPscore, output from cellProximityEnrichment()
#' @param scale scale cell-cell proximity interaction scores
#' @param order_cell_types order cell types based on enrichment correlation
#' @param color_breaks numerical vector of length 3 to represent min, mean
#' and maximum
#' @param color_names character color vector of length 3
#' @returns ggplot heatmap
#' @details This function creates a heatmap that shows the  spatial proximity
#'  enrichment or depletion of cell type pairs.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' x <- cellProximityEnrichment(g, cluster_column = "leiden_clus")
#'
#' cellProximityHeatmap(gobject = g, CPscore = x)
#' @export
cellProximityHeatmap <- function(
        gobject,
        CPscore,
        scale = TRUE,
        order_cell_types = TRUE,
        color_breaks = NULL,
        color_names = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "cellProximityHeatmap") {
    enrich_res <- CPscore$enrichm_res

    # data.table variables
    first_type <- second_type <- unified_int <- NULL

    enrich_res[, first_type := strsplit(
        x = as.character(unified_int), split = "--"
    )[[1]][1],
    by = seq_len(nrow(enrich_res))
    ]
    enrich_res[, second_type := strsplit(
        x = as.character(unified_int), split = "--"
    )[[1]][2],
    by = seq_len(nrow(enrich_res))
    ]

    # create matrix
    enrich_mat <- data.table::dcast.data.table(
        data = enrich_res,
        formula = first_type ~ second_type,
        value.var = "enrichm"
    )
    matrix_d <- as.matrix(enrich_mat[, -1])
    rownames(matrix_d) <- as.vector(enrich_mat[[1]])
    t_matrix_d <- t_flex(matrix_d)

    # fill in NAs based on values in upper and lower matrix triangle
    t_matrix_d[upper.tri(t_matrix_d)][is.na(t_matrix_d[
        upper.tri(t_matrix_d)
    ])] <- matrix_d[upper.tri(matrix_d)][
        is.na(t_matrix_d[upper.tri(t_matrix_d)])
    ]
    t_matrix_d[lower.tri(t_matrix_d)][is.na(t_matrix_d[
        lower.tri(t_matrix_d)
    ])] <- matrix_d[lower.tri(matrix_d)][
        is.na(t_matrix_d[lower.tri(t_matrix_d)])
    ]
    t_matrix_d[is.na(t_matrix_d)] <- 0
    final_matrix <- t_matrix_d

    # scale data
    if (scale == TRUE) {
        final_matrix <- t_flex(scale(t_flex(final_matrix)))
        final_matrix <- t_flex(final_matrix)
        final_matrix[lower.tri(final_matrix)] <- t_flex(final_matrix)[
            lower.tri(final_matrix)
        ]
    }

    # order cell types
    if (order_cell_types == TRUE) {
        cordist <- stats::as.dist(1 - cor_flex(final_matrix))
        clus <- stats::hclust(cordist)
        myorder <- clus$order
        mylabels <- clus$labels
        names(mylabels) <- seq_along(mylabels)
        sample_order <- mylabels[myorder]

        final_matrix <- final_matrix[sample_order, sample_order]
    }

    # create custom colors or not
    if (!is.null(color_breaks) & !is.null(color_names)) {
        if (length(color_breaks) != 3 | !is.numeric(color_breaks)) {
            stop("color_breaks needs to be a numerical vector of length 3")
        }

        if (length(color_names) != 3 | !is.character(color_names)) {
            stop("color_names needs to be a character vector of length 3")
        }

        heatm <- ComplexHeatmap::Heatmap(
            matrix = final_matrix,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            col = GiottoVisuals::colorRamp2(
                breaks = color_breaks, colors = color_names
            )
        )
    } else {
        heatm <- ComplexHeatmap::Heatmap(
            matrix = final_matrix,
            cluster_rows = FALSE,
            cluster_columns = FALSE
        )
    }

    return(plot_output_handler(
        gobject = gobject,
        plot_object = heatm,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}


#' @title cellProximityNetwork
#' @name cellProximityNetwork
#' @description Create network from cell-cell proximity scores
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param CPscore CPscore, output from cellProximityEnrichment()
#' @param remove_self_edges remove enrichment/depletion edges with itself
#' @param self_loop_strength size of self-loops
#' @param color_depletion color for depleted cell-cell interactions
#' @param color_enrichment color for enriched cell-cell interactions
#' @param rescale_edge_weights rescale edge weights (boolean)
#' @param edge_weight_range_depletion numerical vector of length 2 to rescale
#' depleted edge weights
#' @param edge_weight_range_enrichment numerical vector of length 2 to rescale
#' enriched edge weights
#' @param layout layout algorithm to use to draw nodes and edges
#' @param only_show_enrichment_edges show only the enriched pairwise scores
#' @param edge_width_range range of edge width
#' @param node_size size of nodes
#' @param node_color_code color code for the nodes (e.g. cell labels)
#' @param node_text_size size of node labels
#' @returns igraph plot
#' @details This function creates a network that shows the  spatial proximity
#'  enrichment or depletion of cell type pairs.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' x <- cellProximityEnrichment(g, cluster_column = "leiden_clus")
#'
#' cellProximityNetwork(gobject = g, CPscore = x)
#' @export
cellProximityNetwork <- function(
        gobject,
        CPscore,
        remove_self_edges = FALSE,
        self_loop_strength = 0.1,
        color_depletion = "lightgreen",
        color_enrichment = "red",
        rescale_edge_weights = TRUE,
        edge_weight_range_depletion = c(0.1, 1),
        edge_weight_range_enrichment = c(1, 5),
        layout = c("Fruchterman", "DrL", "Kamada-Kawai"),
        only_show_enrichment_edges = FALSE,
        edge_width_range = c(0.1, 2),
        node_size = 4,
        node_color_code = NULL,
        node_text_size = 6,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "cellProximityNetwork") {
    # extract scores

    # data.table variables
    cell_1 <- cell_2 <- unified_int <- color <- size <- name <- NULL

    CPscores <- CPscore[["enrichm_res"]]
    CPscores[, cell_1 := strsplit(
        as.character(unified_int),
        split = "--"
    )[[1]][1],
    by = seq_len(nrow(CPscores))
    ]
    CPscores[, cell_2 := strsplit(
        as.character(unified_int),
        split = "--"
    )[[1]][2],
    by = seq_len(nrow(CPscores))
    ]

    # create igraph with enrichm as weight edges
    igd <- igraph::graph_from_data_frame(
        d = CPscores[, c("cell_1", "cell_2", "enrichm")], directed = FALSE
    )

    if (remove_self_edges == TRUE) {
        igd <- igraph::simplify(
            graph = igd, remove.loops = TRUE, remove.multiple = FALSE
        )
    }

    edges_sizes <- igraph::get.edge.attribute(igd, "enrichm")
    post_edges_sizes <- edges_sizes[edges_sizes > 0]
    neg_edges_sizes <- edges_sizes[edges_sizes <= 0]

    # rescale if wanted
    if (rescale_edge_weights == TRUE) {
        pos_edges_sizes_resc <- scales::rescale(
            x = post_edges_sizes, to = edge_weight_range_enrichment
        )
        neg_edges_sizes_resc <- scales::rescale(
            x = neg_edges_sizes, to = edge_weight_range_depletion
        )
        edges_sizes_resc <- c(pos_edges_sizes_resc, neg_edges_sizes_resc)
    } else {
        edges_sizes_resc <- c(post_edges_sizes, neg_edges_sizes)
    }

    # colors
    edges_colors <- ifelse(edges_sizes > 0, "enriched", "depleted")


    # create coordinates for layout
    if (inherits(layout, c("data.frame", "data.table"))) {
        if (ncol(layout) < 2) {
            stop("custom layout needs to have at least 2 columns")
        }

        if (nrow(layout) != length(igraph::E(igd))) {
            stop("rows of custom layout need to be the same as number of edges")
        }
    } else {
        layout <- match.arg(
            arg = layout, choices = c("Fruchterman", "DrL", "Kamada-Kawai")
        )
    }

    igd <- igraph::set.edge.attribute(
        graph = igd, index = igraph::E(igd), name = "color",
        value = edges_colors
    )
    igd <- igraph::set.edge.attribute(
        graph = igd, index = igraph::E(igd), name = "size",
        value = as.numeric(edges_sizes_resc)
    )

    ## only show attractive edges
    if (only_show_enrichment_edges == TRUE) {
        colors <- igraph::get.edge.attribute(igd, name = "color")
        subvertices_ids <- which(colors == "enriched")
        igd <- igraph::subgraph.edges(graph = igd, eids = subvertices_ids)

        # get new rescale vector (in case vector id is lost)
        edges_sizes_resc <- igraph::E(igd)$size
    }

    ## get coordinates layouts
    if (layout == "Fruchterman") {
        coords <- igraph::layout_with_fr(
            graph = igd, weights = edges_sizes_resc
        )
    } else if (layout == "DrL") {
        coords <- igraph::layout_with_drl(
            graph = igd, weights = edges_sizes_resc
        )
    } else if (layout == "Kamada-Kawai") {
        coords <- igraph::layout_with_kk(
            graph = igd, weights = edges_sizes_resc
        )
    } else {
        stop("Currently no other layouts have been implemented")
    }

    ## create plot
    gpl <- ggraph::ggraph(graph = igd, layout = coords)
    gpl <- gpl + ggraph::geom_edge_link(
        ggplot2::aes(
            color = factor(color),
            edge_width = size, edge_alpha = size
        ),
        show.legend = FALSE
    )

    if (remove_self_edges == FALSE) {
        gpl <- gpl + ggraph::geom_edge_loop(
            ggplot2::aes(
                color = factor(color), edge_width = size,
                edge_alpha = size, strength = self_loop_strength
            ),
            show.legend = FALSE
        )
    }

    gpl <- gpl + ggraph::scale_edge_color_manual(
        values = c("enriched" = color_enrichment, "depleted" = color_depletion)
    )
    gpl <- gpl + ggraph::scale_edge_width(range = edge_width_range)
    gpl <- gpl + ggraph::scale_edge_alpha(range = c(0.1, 1))
    gpl <- gpl + ggraph::geom_node_text(
        ggplot2::aes(label = name),
        repel = TRUE, size = node_text_size
    )
    gpl <- gpl + ggraph::geom_node_point(
        ggplot2::aes(color = name),
        size = node_size
    )
    if (!is.null(node_color_code)) {
        gpl <- gpl + ggplot2::scale_color_manual(values = node_color_code)
    }
    gpl <- gpl + ggplot2::theme_bw() + ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
    )

    return(plot_output_handler(
        gobject = gobject,
        plot_object = gpl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}





#' @title cellProximityVisPlot internals
#' @name cellProximityVisPlot_internals
#' @description
#' Create the plots for `cellProximityVisPlot()`
#' @seealso [cellProximityVisPlot()] [cellProximitySpatPlot3D()]
#' @returns cell proximity plot
NULL


#' @describeIn cellProximityVisPlot_internals Visualize 2D cell-cell
#' interactions according to spatial coordinates in ggplot mode
#' @keywords internal
.cellProximityVisPlot_2D_ggplot <- function(
        gobject,
        interaction_name = NULL,
        cluster_column = NULL,
        sdimx = NULL,
        sdimy = NULL,
        cell_color = NULL,
        cell_color_code = NULL,
        color_as_factor = TRUE,
        show_other_cells = FALSE,
        show_network = FALSE,
        show_other_network = FALSE,
        network_color = NULL,
        spatial_network_name = "Delaunay_network",
        show_grid = FALSE,
        grid_color = NULL,
        spatial_grid_name = "spatial_grid",
        coord_fix_ratio = 1,
        show_legend = TRUE,
        point_size_select = 2,
        point_select_border_col = "black",
        point_select_border_stroke = 0.05,
        point_size_other = 1,
        point_alpha_other = 0.3,
        point_other_border_col = "lightgrey",
        point_other_border_stroke = 0.01,
        ...) {
    # data.table variables
    unified_int <- sdimx_begin <- sdimy_begin <- sdimx_end <- sdimy_end <-
        x_start <- x_end <- NULL
    y_start <- y_end <- cell_ID <- NULL

    if (is.null(interaction_name)) {
        stop("you need to specific at least one interaction name, run
            cellProximityEnrichment")
    }


    cell_locations <- getSpatialLocations(
        gobject = gobject,
        output = "data.table"
    )
    spatial_grid <- getSpatialGrid(
        gobject = gobject,
        name = spatial_grid_name
    )
    cell_metadata <- getCellMetadata(
        gobject = gobject,
        output = "data.table"
    )



    spatial_network <- annotateSpatialNetwork(
        gobject = gobject,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column
    )

    cell_IDs_to_keep <- unique(c(
        spatial_network[unified_int %in% interaction_name]$to,
        spatial_network[unified_int %in% interaction_name]$from
    ))


    if (show_other_cells) {
        CellType <- strsplit(interaction_name, "--")
        all_cell_IDs <- cell_metadata[cell_metadata[[
            cluster_column
        ]] == CellType[[1]][1] |
            cell_metadata[[cluster_column]] == CellType[[1]][2], ]$cell_ID
        other_cell_IDs <- setdiff(all_cell_IDs, cell_IDs_to_keep)
    }


    # annotated cell data
    if (nrow(cell_metadata) == 0) {
        cell_locations_metadata <- cell_locations
    } else {
        cell_locations_metadata <- merge(
            cell_locations, cell_metadata,
            by = "cell_ID"
        )
    }


    # first 2 dimensions need to be defined
    if (is.null(sdimx) | is.null(sdimy)) {
        message("first and second dimenion need to be defined, default is
                first 2\n")
        sdimx <- "sdimx"
        sdimy <- "sdimy"
    }

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()

    if (!is.null(spatial_network) & show_network == TRUE) {
        if (is.null(network_color)) network_color <- "red"
        if (show_other_network) {
            pl <- pl + ggplot2::geom_segment(
                data = spatial_network[!unified_int %in% interaction_name],
                aes(
                    x = sdimx_begin, y = sdimy_begin, xend = sdimx_end,
                    yend = sdimy_end
                ),
                color = "lightgrey", size = 0.5, alpha = 0.5
            )
        }
        pl <- pl + ggplot2::geom_segment(
            data = spatial_network[unified_int %in% interaction_name],
            aes(
                x = sdimx_begin, y = sdimy_begin, xend = sdimx_end,
                yend = sdimy_end
            ),
            color = network_color, size = 0.5, alpha = 0.5
        )
    }

    if (!is.null(spatial_grid) & show_grid == TRUE) {
        if (is.null(grid_color)) grid_color <- "black"
        pl <- pl + ggplot2::geom_rect(
            data = spatial_grid,
            aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end),
            color = grid_color, fill = NA
        )
    }

    # cell color default
    if (is.null(cell_color)) {
        cell_color <- "lightblue"
        pl <- pl + ggplot2::geom_point(
            data = cell_locations[!cell_ID %in% cell_IDs_to_keep],
            aes_string(x = sdimx, y = sdimy),
            show.legend = show_legend, shape = 21, fill = "lightgrey",
            size = point_size_other
        )
        pl <- pl + ggplot2::geom_point(
            data = cell_locations[cell_ID %in% cell_IDs_to_keep],
            aes_string(x = sdimx, y = sdimy),
            show.legend = show_legend, shape = 21, fill = cell_color,
            size = point_size_select
        )
        if (show_other_cells) {
            pl <- pl + ggplot2::geom_point(
                data = cell_locations[cell_ID %in% other_cell_IDs],
                aes_string(x = sdimx, y = sdimy),
                show.legend = show_legend, shape = 21, fill = cell_color,
                alpha = point_alpha_other,
                size = point_size_select * 0.5
            )
        }
    } else if (is.character(cell_color)) {
        if (cell_color %in% colnames(cell_locations_metadata)) {
            if (color_as_factor == TRUE) {
                factor_data <- factor(cell_locations_metadata[[cell_color]])
                cell_locations_metadata[[cell_color]] <- factor_data
            }

            pl <- pl + ggplot2::geom_point(
                data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
                aes_string(x = sdimx, y = sdimy),
                fill = "lightgrey", shape = 21, size = point_size_other,
                color = point_other_border_col,
                stroke = point_other_border_stroke
            )
            pl <- pl + ggplot2::geom_point(
                data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                aes_string(x = sdimx, y = sdimy, fill = cell_color),
                show.legend = show_legend, shape = 21, size = point_size_select,
                color = point_select_border_col,
                stroke = point_select_border_stroke
            )
            if (show_other_cells) {
                pl <- pl + ggplot2::geom_point(
                    data = cell_locations_metadata[cell_ID %in% other_cell_IDs],
                    aes_string(x = sdimx, y = sdimy, fill = cell_color),
                    show.legend = show_legend, shape = 21,
                    alpha = point_alpha_other,
                    size = point_size_select * 0.5
                )
            }



            if (!is.null(cell_color_code)) {
                pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
            } else if (color_as_factor == TRUE) {
                number_colors <- length(unique(factor_data))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instructions(gobject)
                )(n = number_colors)
                names(cell_color_code) <- unique(factor_data)
                pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
            } else if (color_as_factor == FALSE) {
                pl <- pl + set_default_color_continuous_cell(
                    colors = cell_color_code,
                    instrs = instructions(gobject),
                    midpoint = NULL,
                    style = "sequential",
                    type_default = list(
                        pal = c("blue", "red")
                    )
                )
            }
        } else {
            pl <- pl + ggplot2::geom_point(
                data = cell_locations_metadata[
                    !cell_ID %in% cell_IDs_to_keep
                ],
                aes_string(x = sdimx, y = sdimy),
                show.legend = show_legend, shape = 21,
                fill = "lightgrey", size = point_size_other,
                color = point_other_border_col,
                stroke = point_other_border_stroke
            )
            pl <- pl + ggplot2::geom_point(
                data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                aes_string(x = sdimx, y = sdimy),
                show.legend = show_legend, shape = 21, fill = cell_color,
                size = point_size_select,
                color = point_select_border_col,
                stroke = point_select_border_stroke
            )
        }
    }

    pl <- pl + ggplot2::theme_bw() + ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)
    )

    # fix coord ratio
    if (!is.null(coord_fix_ratio)) {
        pl <- pl + ggplot2::coord_fixed(ratio = coord_fix_ratio)
    }

    pl <- pl + ggplot2::labs(x = "x coordinates", y = "y coordinates")

    return(pl)
}



#' @describeIn cellProximityVisPlot_internals Visualize 2D cell-cell
#' interactions according to spatial coordinates in plotly mode
#' @keywords internal
.cellProximityVisPlot_2D_plotly <- function(
        gobject,
        interaction_name = NULL,
        cluster_column = NULL,
        sdimx = NULL,
        sdimy = NULL,
        cell_color = NULL,
        cell_color_code = NULL,
        color_as_factor = TRUE,
        show_other_cells = FALSE,
        show_network = FALSE,
        show_other_network = FALSE,
        network_color = NULL,
        spatial_network_name = "Delaunay_network",
        show_grid = FALSE,
        grid_color = NULL,
        spatial_grid_name = "spatial_grid",
        show_legend = TRUE,
        point_size_select = 2,
        point_size_other = 1,
        point_alpha_other = 0.3,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        ...) {
    package_check("plotly")

    # data.table variables
    cell_ID <- unified_int <- NULL

    if (is.null(interaction_name)) {
        stop("you need to specific at least one interaction name, run
            cellProximityEnrichment")
    }


    cell_locations <- gobject@spatial_locs
    spatial_grid <- gobject@spatial_grid[[spatial_grid_name]]
    cell_metadata <- gobject@cell_metadata


    spatial_network <- annotateSpatialNetwork(
        gobject = gobject,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column
    )

    cell_IDs_to_keep <- unique(c(
        spatial_network[unified_int %in% interaction_name]$to,
        spatial_network[unified_int %in% interaction_name]$from
    ))

    if (show_other_cells) {
        CellType <- strsplit(interaction_name, "-")
        all_cell_IDs <- cell_metadata[cell_metadata[[
            cluster_column
        ]] == CellType[[1]][1] |
            cell_metadata[[cluster_column]] == CellType[[1]][2], ]$cell_ID
        other_cell_IDs <- setdiff(all_cell_IDs, cell_IDs_to_keep)
    }

    # annotated cell data
    if (nrow(cell_metadata) == 0) {
        cell_locations_metadata <- cell_locations
    } else {
        cell_locations_metadata <- merge(
            cell_locations, cell_metadata,
            by = "cell_ID"
        )
    }




    # first 2 dimensions need to be defined
    if (is.null(sdimx) | is.null(sdimy)) {
        message("first and second dimenion need to be defined, default is
                first 2")
        sdimx <- "sdimx"
        sdimy <- "sdimy"
    }


    # create 2D plotly plot

    axis_scale <- match.arg(axis_scale, c("cube", "real", "custom"))

    ratio <- plotly_axis_scale_2D(cell_locations_metadata,
        sdimx = sdimx, sdimy = sdimy,
        mode = axis_scale, custom_ratio = custom_ratio
    )

    pl <- plotly::plot_ly()

    if (show_network == TRUE) {
        if (is.null(network_color)) {
            network_color <- "red"
        }
        if (show_other_network) {
            pl <- pl %>% plotly::add_segments(
                name = paste("unselected", spatial_network_name, sep = " "),
                type = "scatter",
                data = spatial_network[!unified_int %in% interaction_name],
                x = ~sdimx_begin,
                y = ~sdimy_begin,
                xend = ~sdimx_end,
                yend = ~sdimy_end,
                line = list(
                    color = "lightgrey",
                    width = 0.5
                ),
                opacity = 0.3
            )
        }
        pl <- pl %>% plotly::add_segments(
            name = spatial_network_name,
            type = "scatter",
            data = spatial_network[unified_int %in% interaction_name],
            x = ~sdimx_begin,
            y = ~sdimy_begin,
            xend = ~sdimx_end,
            yend = ~sdimy_end,
            line = list(
                color = network_color,
                width = 0.5
            ),
            opacity = 0.8
        )
    }
    if (show_grid == TRUE) {
        if (is.null(grid_color)) {
            grid_color <- "black"
        }
        edges <- plotly_grid(spatial_grid)
        spl <- spl %>% plotly::add_segments(
            name = "spatial_grid",
            type = "scatter",
            data = edges,
            x = ~x,
            y = ~y,
            xend = ~x_end,
            yend = ~y_end,
            line = list(
                color = grid_color,
                width = 1
            ),
            opacity = 1
        )
    }


    if (!is.null(cell_color)) {
        if (cell_color %in% colnames(cell_locations_metadata)) {
            if (is.null(cell_color_code)) {
                number_colors <- length(unique(cell_locations_metadata[[
                    cell_color
                ]]))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instructions(gobject)
                )(n = number_colors)
            }
            cell_locations_metadata[[cell_color]] <- as.factor(
                cell_locations_metadata[[cell_color]]
            )

            pl <- pl %>% plotly::add_trace(
                type = "scatter", mode = "markers",
                # name = "selected cells",
                data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                x = ~sdimx, y = ~sdimy,
                color = cell_locations_metadata[
                    cell_ID %in% cell_IDs_to_keep
                ][[cell_color]],
                colors = cell_color_code,
                marker = list(size = point_size_select)
            )
            if (show_other_cells) {
                pl <- pl %>% plotly::add_trace(
                    type = "scatter", mode = "markers",
                    # name = "selected cells outside network",
                    data = cell_locations_metadata[cell_ID %in% other_cell_IDs],
                    x = ~sdimx, y = ~sdimy,
                    color = cell_locations_metadata[
                        cell_ID %in% other_cell_IDs
                    ][[cell_color]],
                    colors = cell_color_code,
                    opacity = point_alpha_other,
                    marker = list(size = point_size_select * 0.7)
                )
            }
            pl <- pl %>% plotly::add_trace(
                type = "scatter", mode = "markers",
                name = "unselected cells",
                data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
                x = ~sdimx, y = ~sdimy,
                marker = list(
                    size = point_size_other,
                    color = "lightgray",
                    colors = "lightgray"
                ),
                opacity = point_alpha_other
            )
        } else {
            message("cell_color not exist!")
        }
    } else {
        pl <- pl %>% plotly::add_trace(
            type = "scatter", mode = "markers",
            name = "selected cells",
            data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
            x = ~sdimx, y = ~sdimy,
            marker = list(
                size = point_size_select,
                color = "lightblue",
                colors = "lightblue"
            )
        )
        if (show_other_cells) {
            pl <- pl %>% plotly::add_trace(
                type = "scatter", mode = "markers",
                data = cell_locations_metadata[cell_ID %in% other_cell_IDs],
                x = ~sdimx, y = ~sdimy,
                name = "selected cells outside network",
                marker = list(
                    size = point_size_select * 0.7,
                    color = "lightblue",
                    colors = "lightblue"
                ),
                opacity = point_alpha_other
            )
        }
        pl <- pl %>% plotly::add_trace(
            type = "scatter", mode = "markers",
            name = "unselected cells",
            data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
            x = ~sdimx, y = ~sdimy,
            marker = list(
                size = point_size_other,
                color = "lightgray",
                colors = "lightgray"
            ),
            opacity = point_alpha_other
        )
    }

    pl <- pl %>%
        plotly::layout(
            list(
                xaxis = list(title = "X", nticks = x_ticks),
                yaxis = list(title = "Y", nticks = y_ticks)
            ),
            legend = list(
                x = 100, y = 0.5,
                font = list(family = "sans-serif", size = 12)
            )
        )
    return((pl))
}


#' @describeIn cellProximityVisPlot_internals Visualize 3D cell-cell
#' interactions according to spatial coordinates in plotly mode
#' @keywords internal
.cellProximityVisPlot_3D_plotly <- function(
        gobject,
        interaction_name = NULL,
        cluster_column = NULL,
        sdimx = NULL,
        sdimy = NULL,
        sdimz = NULL,
        cell_color = NULL,
        cell_color_code = NULL,
        color_as_factor = TRUE,
        show_other_cells = FALSE,
        show_network = FALSE,
        show_other_network = FALSE,
        network_color = NULL,
        spatial_network_name = "Delaunay_network",
        show_grid = FALSE,
        grid_color = NULL,
        spatial_grid_name = "spatial_grid",
        show_legend = TRUE,
        point_size_select = 2,
        point_size_other = 1,
        point_alpha_other = 0.5,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        z_ticks = NULL,
        ...) {
    package_check("plotly")

    # data.table variables
    cell_ID <- unified_int <- NULL

    if (is.null(interaction_name)) {
        stop("you need to specific at least one interaction name,
            run cellProximityEnrichment")
    }


    cell_locations <- gobject@spatial_locs
    spatial_grid <- gobject@spatial_grid[[spatial_grid_name]]
    cell_metadata <- gobject@cell_metadata


    spatial_network <- annotateSpatialNetwork(
        gobject = gobject,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column
    )

    cell_IDs_to_keep <- unique(c(
        spatial_network[
            unified_int %in% interaction_name
        ]$to,
        spatial_network[unified_int %in% interaction_name]$from
    ))

    if (show_other_cells) {
        CellType <- strsplit(interaction_name, "-")
        all_cell_IDs <- cell_metadata[cell_metadata[[
            cluster_column
        ]] == CellType[[1]][1] |
            cell_metadata[[cluster_column]] == CellType[[1]][2], ]$cell_ID
        other_cell_IDs <- setdiff(all_cell_IDs, cell_IDs_to_keep)
    }

    # annotated cell data
    if (nrow(cell_metadata) == 0) {
        cell_locations_metadata <- cell_locations
    } else {
        cell_locations_metadata <- merge(
            cell_locations, cell_metadata,
            by = "cell_ID"
        )
    }




    # first 2 dimensions need to be defined
    if (is.null(sdimx) | is.null(sdimy)) {
        message("first and second dimenion need to be defined, default is
                first 2")
        sdimx <- "sdimx"
        sdimy <- "sdimy"
    }





    # if 3 dimensions are defined create a 3D plot

    # create 3D plotly plot

    pl <- plotly::plot_ly()

    axis_scale <- match.arg(axis_scale, c("cube", "real", "custom"))

    ratio <- plotly_axis_scale_3D(cell_locations_metadata,
        sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
        mode = axis_scale, custom_ratio = custom_ratio
    )

    if (!is.null(cell_color)) {
        if (cell_color %in% colnames(cell_locations_metadata)) {
            if (is.null(cell_color_code)) {
                number_colors <- length(unique(cell_locations_metadata[[
                    cell_color
                ]]))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instructions(gobject)
                )(n = number_colors)
            }
            cell_locations_metadata[[cell_color]] <- as.factor(
                cell_locations_metadata[[cell_color]]
            )

            pl <- pl %>%
                plotly::add_trace(
                    type = "scatter3d", mode = "markers",
                    data = cell_locations_metadata[
                        cell_ID %in% cell_IDs_to_keep
                    ],
                    x = ~sdimx, y = ~sdimy, z = ~sdimz,
                    color = cell_locations_metadata[
                        cell_ID %in% cell_IDs_to_keep
                    ][[cell_color]],
                    colors = cell_color_code,
                    marker = list(size = point_size_select)
                ) %>%
                plotly::add_trace(
                    type = "scatter3d", mode = "markers",
                    name = "unselected cells",
                    data = cell_locations_metadata[
                        !cell_ID %in% cell_IDs_to_keep
                    ],
                    x = ~sdimx, y = ~sdimy, z = ~sdimz,
                    marker = list(
                        size = point_size_other,
                        color = "lightgray",
                        colors = "lightgray"
                    ),
                    opacity = point_alpha_other
                )
            if (show_other_cells) {
                pl <- pl %>% plotly::add_trace(
                    type = "scatter3d", mode = "markers",
                    # name = "selected cells outside network",
                    data = cell_locations_metadata[cell_ID %in% other_cell_IDs],
                    x = ~sdimx, y = ~sdimy, z = ~sdimz,
                    color = cell_locations_metadata[
                        cell_ID %in% other_cell_IDs
                    ][[cell_color]],
                    colors = cell_color_code,
                    opacity = point_alpha_other,
                    marker = list(size = point_size_select * 0.7)
                )
            }
        } else {
            message("cell_color not exist!")
        }
    } else {
        pl <- pl %>%
            plotly::add_trace(
                type = "scatter3d", mode = "markers",
                name = "selected cells",
                data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                x = ~sdimx, y = ~sdimy, z = ~sdimz,
                marker = list(
                    size = point_size_select,
                    color = "lightblue",
                    colors = "lightblue"
                )
            ) %>%
            plotly::add_trace(
                type = "scatter3d", mode = "markers",
                name = "unselected cells",
                data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
                x = ~sdimx, y = ~sdimy, z = ~sdimz,
                marker = list(
                    size = point_size_other,
                    color = "lightgray",
                    colors = "lightgray"
                ),
                opacity = point_alpha_other
            )
        if (show_other_cells) {
            pl <- pl %>% plotly::add_trace(
                type = "scatter3d", mode = "markers",
                data = cell_locations_metadata[cell_ID %in% other_cell_IDs],
                x = ~sdimx, y = ~sdimy, z = ~sdimz,
                name = "selected cells outside network",
                marker = list(
                    size = point_size_select * 0.7,
                    color = "lightblue",
                    colors = "lightblue"
                ),
                opacity = point_alpha_other
            )
        }
    }
    if (!is.null(spatial_network) & show_network == TRUE) {
        if (is.null(network_color)) {
            network_color <- "red"
        }
        unselect_network <- spatial_network[!unified_int %in% interaction_name]
        select_network <- spatial_network[unified_int %in% interaction_name]
        pl <- pl %>% plotly::add_trace(
            name = "sptial network", mode = "lines",
            type = "scatter3d", opacity = 0.5,
            data = plotly_network(select_network),
            x = ~x, y = ~y, z = ~z, inherit = FALSE,
            line = list(color = network_color)
        )
        if (show_other_network == TRUE) {
            pl <- pl %>% plotly::add_trace(
                name = "unselected sptial network", mode = "lines",
                type = "scatter3d", opacity = 0.1,
                data = plotly_network(unselect_network),
                x = ~x, y = ~y, z = ~z, inherit = FALSE,
                line = list(color = "lightgray")
            )
        }
    }

    pl <- pl %>% plotly::layout(scene = list(
        xaxis = list(title = "X", nticks = x_ticks),
        yaxis = list(title = "Y", nticks = y_ticks),
        zaxis = list(title = "Z", nticks = z_ticks),
        aspectmode = "manual",
        aspectratio = list(
            x = ratio[[1]],
            y = ratio[[2]],
            z = ratio[[3]]
        )
    ))
    return(pl)
}

#' @title cellProximityVisPlot
#' @name cellProximityVisPlot
#' @description Visualize cell-cell interactions according to spatial
#' coordinates
#' @param gobject giotto object
#' @param interaction_name cell-cell interaction name
#' @param cluster_column cluster column with cell clusters
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimz')
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param show_other_cells show not selected cells
#' @param show_network show underlying spatial network
#' @param show_other_network show underlying spatial network of other cells
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param show_legend show legend
#' @param point_size_select size of selected points
#' @param point_select_border_col border color of selected points
#' @param point_select_border_stroke stroke size of selected points
#' @param point_size_other size of other points
#' @param point_alpha_other alpha of other points
#' @param point_other_border_col border color of other points
#' @param point_other_border_stroke stroke size of other points
#' @param axis_scale scale of axis
#' @param custom_ratio custom ratio of scales
#' @param x_ticks x ticks
#' @param y_ticks y ticks
#' @param z_ticks z ticks
#' @param plot_method method to plot
#' @param \dots additional parameters
#' @returns ggplot or plotly
#' @details Description of parameters.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- createSpatialGrid(g, sdimx_stepsize = 5, sdimy_stepsize = 5)
#' x <- cellProximityEnrichment(g, cluster_column = "leiden_clus")
#'
#' cellProximityVisPlot(
#'     gobject = g, interaction_name = x,
#'     cluster_column = "leiden_clus", sdimx = "sdimx", sdimy = "sdimy"
#' )
#' @export
cellProximityVisPlot <- function(
        gobject,
        interaction_name = NULL,
        cluster_column = NULL,
        sdimx = NULL,
        sdimy = NULL,
        sdimz = NULL,
        cell_color = NULL,
        cell_color_code = NULL,
        color_as_factor = TRUE,
        show_other_cells = FALSE,
        show_network = FALSE,
        show_other_network = FALSE,
        network_color = NULL,
        spatial_network_name = "Delaunay_network",
        show_grid = FALSE,
        grid_color = NULL,
        spatial_grid_name = "spatial_grid",
        coord_fix_ratio = 1,
        show_legend = TRUE,
        point_size_select = 2,
        point_select_border_col = "black",
        point_select_border_stroke = 0.05,
        point_size_other = 1,
        point_alpha_other = 0.3,
        point_other_border_col = "lightgrey",
        point_other_border_stroke = 0.01,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        z_ticks = NULL,
        plot_method = c("ggplot", "plotly"),
        ...) {
    ## decide plot method
    plot_method <- match.arg(plot_method, choices = c("ggplot", "plotly"))
    axis_scale <- match.arg(axis_scale, c("cube", "real", "custom"))


    if (plot_method == "ggplot") {
        if (is.null(sdimx) | is.null(sdimy)) {
            warning("plot_method = ggplot, but spatial dimensions for sdimx
            and sdimy for 2D plotting are not given. \n
            It will default to the 'sdimx' and 'sdimy'")
            sdimx <- "sdimx"
            sdimy <- "sdimy"
        }

        if (length(c(sdimx, sdimy, sdimz)) == 3) {
            warning("ggplot is not able to produce 3D plot! Please choose
                    plotly method")
        }
        result <- .cellProximityVisPlot_2D_ggplot(
            gobject = gobject,
            interaction_name = interaction_name,
            cluster_column = cluster_column,
            sdimx = sdimx,
            sdimy = sdimy,
            cell_color = cell_color,
            cell_color_code = cell_color_code,
            color_as_factor = color_as_factor,
            show_other_cells = show_other_cells,
            show_network = show_network,
            show_other_network = show_other_network,
            network_color = network_color,
            spatial_network_name = spatial_network_name,
            show_grid = show_grid,
            grid_color = grid_color,
            spatial_grid_name = spatial_grid_name,
            coord_fix_ratio = coord_fix_ratio,
            show_legend = show_legend,
            point_size_select = point_size_select,
            point_select_border_col = point_select_border_col,
            point_select_border_stroke = point_select_border_stroke,
            point_size_other = point_size_other,
            point_alpha_other = point_alpha_other,
            point_other_border_col = point_other_border_col,
            point_other_border_stroke = point_other_border_stroke,
            ...
        )
    } else if (plot_method == "plotly") {
        if (length(c(sdimx, sdimy, sdimz)) == 3) {
            result <- .cellProximityVisPlot_3D_plotly(
                gobject = gobject,
                interaction_name = interaction_name,
                cluster_column = cluster_column,
                sdimx = sdimx,
                sdimy = sdimy,
                sdimz = sdimz,
                cell_color = cell_color,
                cell_color_code = cell_color_code,
                color_as_factor = color_as_factor,
                show_other_cells = show_other_cells,
                show_network = show_network,
                show_other_network = show_other_network,
                network_color = network_color,
                spatial_network_name = spatial_network_name,
                show_grid = show_grid,
                grid_color = grid_color,
                spatial_grid_name = spatial_grid_name,
                show_legend = show_legend,
                point_size_select = point_size_select,
                point_size_other = point_size_other,
                point_alpha_other = point_alpha_other,
                axis_scale = axis_scale,
                custom_ratio = custom_ratio,
                x_ticks = x_ticks,
                y_ticks = y_ticks,
                z_ticks = z_ticks,
                ...
            )
        } else {
            if (is.null(sdimx) | is.null(sdimy)) {
                stop("plotly in 2D requires you to define sdimx and sdimy")
            }

            ## run: visPlot_2D_plotly
            result <- .cellProximityVisPlot_2D_plotly(
                gobject = gobject,
                interaction_name = interaction_name,
                cluster_column = cluster_column,
                sdimx = sdimx,
                sdimy = sdimy,
                cell_color = cell_color,
                cell_color_code = cell_color_code,
                color_as_factor = color_as_factor,
                show_other_cells = show_other_cells,
                show_network = show_network,
                show_other_network = show_other_network,
                network_color = network_color,
                spatial_network_name = spatial_network_name,
                show_grid = show_grid,
                grid_color = grid_color,
                spatial_grid_name = spatial_grid_name,
                show_legend = show_legend,
                point_size_select = point_size_select,
                point_size_other = point_size_other,
                point_alpha_other = point_alpha_other,
                axis_scale = axis_scale,
                custom_ratio = custom_ratio,
                x_ticks = x_ticks,
                y_ticks = y_ticks,
                ...
            )
        }
    }
    return(result)
}







# * ####
# NEW plots ####




#' @title plotCellProximityFeats
#' @name plotCellProximityFeats
#' @description Create visualization for cell proximity feature scores
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @param icfObject ICF (interaction changed feature) score object
#' @param method plotting method to use
#' @param min_cells minimum number of source cell type
#' @param min_cells_expr minimum expression level for source cell type
#' @param min_int_cells minimum number of interacting neighbor cell type
#' @param min_int_cells_expr minimum expression level for interacting neighbor
#' cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum log2 fold-change
#' @param min_zscore minimum z-score change
#' @param zscores_column calculate z-scores over cell types or features
#' @param direction differential expression directions to keep
#' @returns volcano, cell_barplot, cell-cell, cell_sankey, heatmap, or dotplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' icfObject <- findInteractionChangedFeats(g, cluster_column = "leiden_clus")
#'
#' plotCellProximityFeats(
#'     gobject = g, icfObject = icfObject,
#'     show_plot = TRUE, save_plot = FALSE, return_plot = FALSE
#' )
#' @export
plotCellProximityFeats <- function(
        gobject,
        icfObject,
        method = c(
            "volcano", "cell_barplot", "cell-cell", "cell_sankey",
            "heatmap", "dotplot"
        ),
        min_cells = 4,
        min_cells_expr = 1,
        min_int_cells = 4,
        min_int_cells_expr = 1,
        min_fdr = 0.1,
        min_spat_diff = 0.2,
        min_log2_fc = 0.2,
        min_zscore = 2,
        zscores_column = c("cell_type", "feats"),
        direction = c("both", "up", "down"),
        cell_color_code = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotCellProximityFeats") {
    if (!"icfObject" %in% class(icfObject)) {
        stop("icfObject needs to be the output from
            findInteractionChangedFeats() or findICF()")
    }

    # print, return and save parameters
    show_plot <- ifelse(
        is.na(show_plot),
        readGiottoInstructions(gobject, param = "show_plot"), show_plot
    )
    save_plot <- ifelse(
        is.na(save_plot),
        readGiottoInstructions(gobject, param = "save_plot"), save_plot
    )
    return_plot <- ifelse(
        is.na(return_plot),
        readGiottoInstructions(gobject, param = "return_plot"), return_plot
    )


    ## first filter
    filter_icf <- filterInteractionChangedFeats(
        icfObject = icfObject,
        min_cells = min_cells,
        min_cells_expr = min_cells_expr,
        min_int_cells = min_int_cells,
        min_int_cells_expr = min_int_cells_expr,
        min_fdr = min_fdr,
        min_spat_diff = min_spat_diff,
        min_log2_fc = min_log2_fc,
        min_zscore = min_zscore,
        zscores_column = zscores_column,
        direction = direction
    )

    complete_part <- filter_icf[["ICFscores"]]

    ## other parameters
    method <- match.arg(
        method,
        choices = c(
            "volcano", "cell_barplot", "cell-cell", "cell_sankey",
            "heatmap", "dotplot"
        )
    )


    # variables
    log2fc <- p.adj <- unif_int <- N <- cell_type <- int_cell_type <- NULL

    ## create data.table for visualization
    if (method == "volcano") {
        ## volcanoplot
        pl <- ggplot2::ggplot()
        pl <- pl + ggplot2::geom_point(
            data = complete_part,
            ggplot2::aes(
                x = log2fc,
                y = ifelse(is.infinite(-log10(p.adj)),
                    1000, -log10(p.adj)
                )
            )
        )
        pl <- pl + ggplot2::theme_classic()
        pl <- pl + ggplot2::geom_vline(xintercept = 0, linetype = 2)
        pl <- pl + ggplot2::labs(
            x = "log2 fold-change", y = "-log10(p.adjusted)"
        )


        ## print plot
        if (show_plot == TRUE) {
            print(pl)
        }

        ## save plot
        if (save_plot == TRUE) {
            do.call(
                "all_plots_save_function",
                c(list(
                    gobject = gobject, plot_object = pl,
                    default_save_name = default_save_name
                ), save_param)
            )
        }

        ## return plot
        if (return_plot == TRUE) {
            return(pl)
        }
    } else if (method == "cell-cell") {
        nr_int_selection_scores <- complete_part[, .N, by = unif_int]
        order_interactions <- nr_int_selection_scores[order(N)]$unif_int

        complete_part[, unif_int := factor(unif_int, order_interactions)]

        pl <- ggplot2::ggplot()
        pl <- pl + ggplot2::geom_bar(
            data = complete_part,
            ggplot2::aes(x = unif_int, fill = unif_int)
        )
        pl <- pl + ggplot2::theme_classic() + ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 90, hjust = 1, vjust = 1
            )
        )
        pl <- pl + ggplot2::coord_flip()

        ## print plot
        if (show_plot == TRUE) {
            print(pl)
        }

        ## save plot
        if (save_plot == TRUE) {
            do.call(
                "all_plots_save_function",
                c(list(
                    gobject = gobject, plot_object = pl,
                    default_save_name = default_save_name
                ), save_param)
            )
        }

        ## return plot
        if (return_plot == TRUE) {
            return(pl)
        }
    } else if (method == "cell_barplot") {
        # by source cell type plot
        nr_source_selection_scores <- complete_part[, .N, by = cell_type]
        order_source <- nr_source_selection_scores[order(N)]$cell_type

        complete_part[, cell_type := factor(cell_type, order_source)]

        pl <- ggplot2::ggplot()
        pl <- pl + ggplot2::geom_bar(
            data = complete_part,
            ggplot2::aes(x = cell_type, fill = int_cell_type)
        )
        if (!is.null(cell_color_code)) {
            pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
        }
        pl <- pl + ggplot2::theme_classic() + ggplot2::theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
        )
        pl <- pl + ggplot2::labs(
            x = "", y = "# of features influenced by cell neighborhood"
        )


        ## print plot
        if (show_plot == TRUE) {
            print(pl)
        }

        ## save plot
        if (save_plot == TRUE) {
            do.call(
                "all_plots_save_function",
                c(list(
                    gobject = gobject, plot_object = pl,
                    default_save_name = default_save_name
                ), save_param)
            )
        }

        ## return plot
        if (return_plot == TRUE) {
            return(pl)
        }
    } else if (method == "cell_sankey") {
        testalluv <- complete_part[, .N, by = c("int_cell_type", "cell_type")]

        # package check for ggalluvial
        # verify if optional package is installed
        package_check(pkg_name = "ggalluvial", repository = "CRAN")

        pl <- ggplot2::ggplot(
            testalluv,
            ggplot2::aes(y = N, axis1 = cell_type, axis2 = int_cell_type)
        ) +
            ggalluvial::geom_alluvium(aes(fill = cell_type), width = 1 / 12) +
            ggalluvial::geom_stratum(
                width = 1 / 12, fill = "black", color = "grey"
            ) +
            ggplot2::scale_x_discrete(
                limits = c("cell type", "neighbours"), expand = c(.05, .05)
            ) +
            ggplot2::geom_label(
                stat = "stratum", label.strata = TRUE, size = 3
            ) +
            ggplot2::theme_classic() +
            ggplot2::labs(
                x = "", y = "# of features influenced by cell neighborhood"
            )

        if (!is.null(cell_color_code)) {
            pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
        }



        ## print plot
        if (show_plot == TRUE) {
            print(pl)
        }

        ## save plot
        if (save_plot == TRUE) {
            do.call(
                "all_plots_save_function",
                c(list(
                    gobject = gobject, plot_object = pl,
                    default_save_name = default_save_name
                ), save_param)
            )
        }

        ## return plot
        if (return_plot == TRUE) {
            return(pl)
        }
    } else if (method == "dotplot") {
        changed_feats <- complete_part[
            , .N,
            by = c("cell_type", "int_cell_type")
        ]

        changed_feats[, cell_type := factor(cell_type, unique(cell_type))]
        changed_feats[, int_cell_type := factor(
            int_cell_type, unique(int_cell_type)
        )]

        pl <- ggplot2::ggplot()
        pl <- pl + ggplot2::theme_classic()
        pl <- pl + ggplot2::geom_point(
            data = changed_feats,
            ggplot2::aes(x = cell_type, y = int_cell_type, size = N)
        )
        pl <- pl + ggplot2::scale_size_continuous(
            guide = guide_legend(title = "# of ICFs")
        )
        pl <- pl + ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90, vjust = 1, hjust = 1
        ))
        pl <- pl + ggplot2::labs(
            x = "source cell type", y = "neighbor cell type"
        )

        ## print plot
        if (show_plot == TRUE) {
            print(pl)
        }

        ## save plot
        if (save_plot == TRUE) {
            do.call(
                "all_plots_save_function",
                c(list(
                    gobject = gobject, plot_object = pl,
                    default_save_name = default_save_name
                ), save_param)
            )
        }

        ## return plot
        if (return_plot == TRUE) {
            return(pl)
        }
    } else if (method == "heatmap") {
        changed_feats <- complete_part[
            , .N,
            by = c("cell_type", "int_cell_type")
        ]

        changed_feats[, cell_type := factor(cell_type, unique(cell_type))]
        changed_feats[, int_cell_type := factor(
            int_cell_type, unique(int_cell_type)
        )]

        changed_feats_d <- data.table::dcast.data.table(
            changed_feats, cell_type ~ int_cell_type,
            value.var = "N", fill = 0
        )
        changed_feats_m <- dt_to_matrix(changed_feats_d)

        col_fun <- GiottoVisuals::colorRamp2(
            breaks = stats::quantile(log2(changed_feats_m + 1)),
            colors = c("white", "white", "blue", "yellow", "red")
        )

        heatm <- ComplexHeatmap::Heatmap(log2(changed_feats_m + 1),
            col = col_fun,
            row_title = "cell_type",
            column_title = "int_cell_type",
            heatmap_legend_param = list(title = "log2(# DEGs)")
        )

        ## print plot
        if (show_plot == TRUE) {
            print(heatm)
        }

        ## save plot
        if (save_plot == TRUE) {
            do.call(
                "all_plots_save_function",
                c(list(
                    gobject = gobject, plot_object = heatm,
                    default_save_name = default_save_name
                ), save_param)
            )
        }

        ## return plot
        if (return_plot == TRUE) {
            return(heatm)
        }
    }
}













#' @title plotCPF
#' @name plotCPF
#' @description Create visualization for cell proximity feature scores
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @param icfObject ICF (interaction changed feature) score object
#' @param method plotting method to use
#' @param min_cells minimum number of source cell type
#' @param min_cells_expr minimum expression level for source cell type
#' @param min_int_cells minimum number of interacting neighbor cell type
#' @param min_int_cells_expr minimum expression level for interacting neighbor
#' cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum log2 fold-change
#' @param min_zscore minimum z-score change
#' @param zscores_column calculate z-scores over cell types or features
#' @param direction differential expression directions to keep
#' @returns plot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' icfObject <- findInteractionChangedFeats(g,
#'     cluster_column = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17"), nr_permutations = 10
#' )
#'
#' plotCPF(
#'     gobject = g, icfObject = icfObject, show_plot = TRUE,
#'     save_plot = FALSE, return_plot = FALSE
#' )
#' @export
plotCPF <- function(
        gobject,
        icfObject,
        method = c(
            "volcano", "cell_barplot", "cell-cell", "cell_sankey",
            "heatmap", "dotplot"
        ),
        min_cells = 5,
        min_cells_expr = 1,
        min_int_cells = 3,
        min_int_cells_expr = 1,
        min_fdr = 0.05,
        min_spat_diff = 0.2,
        min_log2_fc = 0.2,
        min_zscore = 2,
        zscores_column = c("cell_type", "feats"),
        direction = c("both", "up", "down"),
        cell_color_code = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotCPG") {
    plotCellProximityFeats(
        gobject = gobject,
        icfObject = icfObject,
        method = method,
        min_cells = min_cells,
        min_cells_expr = min_cells_expr,
        min_int_cells = min_int_cells,
        min_int_cells_expr = min_int_cells_expr,
        min_fdr = min_fdr,
        min_spat_diff = min_spat_diff,
        min_log2_fc = min_log2_fc,
        min_zscore = min_zscore,
        zscores_column = zscores_column,
        direction = direction,
        cell_color_code = cell_color_code,
        show_plot = show_plot,
        return_plot = return_plot,
        save_plot = save_plot,
        save_param = save_param,
        default_save_name = default_save_name
    )
}



#' @title Plot interaction changed features
#' @name plotInteractionChangedFeats
#' @description Create barplot to visualize interaction changed features
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @param icfObject ICF (interaction changed feature) score object
#' @param source_type cell type of the source cell
#' @param source_markers markers for the source cell type
#' @param ICF_feats named character vector of ICF features
#' @returns plot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' icfObject <- findInteractionChangedFeats(g,
#'     cluster_column = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17"), nr_permutations = 10
#' )
#'
#' plotInteractionChangedFeats(
#'     gobject = g, icfObject = icfObject,
#'     source_type = "1", source_markers = "Ccnd2",
#'     ICF_feats = c("3" = "Gna12", "1" = "Ccnd2", "8" = "Btbd17")
#' )
#' @export
plotInteractionChangedFeats <- function(
        gobject,
        icfObject,
        source_type,
        source_markers,
        ICF_feats,
        cell_color_code = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotInteractionChangedFeats") {
    # data.table variables
    cell_type <- int_cell_type <- log2fc <- NULL


    if (!"icfObject" %in% class(icfObject)) {
        stop("icfObject needs to be the output from
            findInteractionChangedFeats() or findICF()")
    }

    ICFscores <- icfObject[["ICFscores"]]

    # combine feats
    names(source_markers) <- rep("marker", length(source_markers))
    neighbor_types <- names(ICF_feats)
    all_feats <- c(source_markers, ICF_feats)

    # warning if there are feats selected that are not detected
    detected_feats <- unique(ICFscores[["feats"]])
    not_detected_feats <- all_feats[!all_feats %in% detected_feats]
    if (length(not_detected_feats) > 0) {
        cat(
            "These selected features are not in the icfObject: \n",
            not_detected_feats
        )
    }

    # data.table set column names
    feats <- group <- NULL

    tempDT <- ICFscores[feats %in% all_feats][cell_type == source_type][
        int_cell_type %in% neighbor_types
    ]
    tempDT[, feats := factor(feats, levels = detected_feats)]
    tempDT[, group := names(ICF_feats[ICF_feats == feats]),
        by = seq_len(nrow(tempDT))
    ]


    if (is.null(cell_color_code)) {
        mycolors <- set_default_color_discrete_cell(
            instrs = instructions(gobject)
        )(n = length(unique(
            tempDT$int_cell_type
        )))
        names(mycolors) <- unique(tempDT$int_cell_type)
    } else {
        mycolors <- cell_color_code
    }


    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(
        axis.text.x = ggplot2::element_text(
            size = 14, angle = 45, vjust = 1, hjust = 1
        ),
        axis.text.y = ggplot2::element_text(size = 14),
        axis.title = ggplot2::element_text(size = 14)
    )
    pl <- pl + ggplot2::geom_bar(
        data = tempDT,
        ggplot2::aes(x = feats, y = log2fc, fill = int_cell_type),
        stat = "identity", position = ggplot2::position_dodge()
    )
    pl <- pl + ggplot2::scale_fill_manual(values = mycolors)
    pl <- pl + ggplot2::labs(x = "", title = paste0(
        "fold-change z-scores in ", source_type
    ))

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}







#' @title Plot interaction changed features
#' @name plotICF
#' @description Create barplot to visualize interaction changed features
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @param icfObject ICF (interaction changed feature) score object
#' @param source_type cell type of the source cell
#' @param source_markers markers for the source cell type
#' @param ICF_feats named character vector of ICF features
#' @returns plot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' icfObject <- findInteractionChangedFeats(g,
#'     cluster_column = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17"), nr_permutations = 10
#' )
#'
#' plotICF(
#'     gobject = g, icfObject = icfObject,
#'     source_type = "1", source_markers = "Ccnd2",
#'     ICF_feats = c("3" = "Gna12", "1" = "Ccnd2", "8" = "Btbd17")
#' )
#' @export
plotICF <- function(
        gobject,
        icfObject,
        source_type,
        source_markers,
        ICF_feats,
        cell_color_code = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotICF") {
    plotInteractionChangedFeats(
        gobject = gobject,
        icfObject = icfObject,
        source_type = source_type,
        source_markers = source_markers,
        ICF_feats = ICF_feats,
        cell_color_code = cell_color_code,
        show_plot = show_plot,
        return_plot = return_plot,
        save_plot = save_plot,
        save_param = save_param,
        default_save_name = default_save_name
    )
}








#' @title plotCombineInteractionChangedFeats
#' @name plotCombineInteractionChangedFeats
#' @description Create visualization for combined (pairwise) ICF scores
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param combIcfObject ICFscores, output from combineInteractionChangedFeats()
#' @param selected_interactions interactions to show
#' @param selected_feat_to_feat pairwise feature combinations to show
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or feats with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param colors vector with two colors to use
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' g_icf <- findInteractionChangedFeats(g,
#'     cluster_column = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17"), nr_permutations = 10
#' )
#'
#' combIcfObject <- combineInteractionChangedFeats(g_icf)
#'
#' plotCombineInteractionChangedFeats(
#'     gobject = g,
#'     combIcfObject = combIcfObject,
#'     selected_feat_to_feat = c("Btbd17--Ccnd2", "Btbd17--Gna12"),
#'     selected_interactions = "1--8"
#' )
#' @export
plotCombineInteractionChangedFeats <- function(
        gobject,
        combIcfObject,
        selected_interactions = NULL,
        selected_feat_to_feat = NULL,
        detail_plot = TRUE,
        simple_plot = FALSE,
        simple_plot_facet = c("interaction", "feats"),
        facet_scales = "fixed",
        facet_ncol = length(selected_feat_to_feat),
        facet_nrow = length(selected_interactions),
        colors = c("#9932CC", "#FF8C00"),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotCombineICF") {
    ## check validity
    if (!"combIcfObject" %in% class(combIcfObject)) {
        stop("combIcfObject needs to be the output from
            combineInteractionChangedFeats() or combineICF()")
    }
    combIcfscore <- copy(combIcfObject[["combICFscores"]])

    if (is.null(selected_interactions) | is.null(selected_feat_to_feat)) {
        stop("You need to provide a selection of cell-cell interactions and
            features-features to plot")
    }


    # data.table variables
    unif_feat_feat <- unif_int <- other_2 <- sel_2 <- other_1 <- sel_1 <-
        cols <- NULL


    subDT <- combIcfscore[unif_feat_feat %in% selected_feat_to_feat &
        unif_int %in% selected_interactions]

    # order interactions and feat-to-feat according to input
    subDT[, unif_feat_feat := factor(
        unif_feat_feat,
        levels = selected_feat_to_feat
    )]
    subDT[, unif_int := factor(unif_int, levels = selected_interactions)]

    if (simple_plot == FALSE) {
        pl <- ggplot2::ggplot()
        pl <- pl + ggplot2::theme_bw()

        if (detail_plot == TRUE) {
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(x = 0, y = other_2, colour = "other cell expression"),
                shape = 1
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(x = 0, y = sel_2, colour = "selected cell expression"),
                shape = 1
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(x = other_1, y = 0, colour = "other cell expression"),
                shape = 1
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(x = sel_1, y = 0, colour = "selected cell expression"),
                shape = 1
            )
        }

        pl <- pl + ggplot2::geom_point(
            data = subDT,
            aes(x = other_1, y = other_2, colour = "other cell expression"),
            size = 2
        )
        pl <- pl + ggplot2::geom_point(
            data = subDT,
            aes(x = sel_1, y = sel_2, colour = "selected cell expression"),
            size = 2
        )
        pl <- pl + ggplot2::geom_segment(data = subDT, aes(
            x = other_1, xend = sel_1,
            y = other_2, yend = sel_2
        ), linetype = 2)
        pl <- pl + ggplot2::labs(
            x = paste(subDT$feats_1, subDT$cell_type_1, sep = " in "),
            y = paste(subDT$feats_2, subDT$cell_type_2, sep = " in ")
        )
        pl <- pl + ggplot2::scale_colour_manual(
            name = "expression source", values = colors
        )
        pl <- pl + ggplot2::facet_wrap(~ unif_feat_feat + unif_int,
            nrow = facet_nrow, ncol = facet_ncol,
            scales = facet_scales
        )
    } else {
        simple_plot_facet <- match.arg(
            arg = simple_plot_facet, choices = c("interaction", "feats")
        )

        if (simple_plot_facet == "interaction") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_bw()
            pl <- pl + ggplot2::geom_segment(data = subDT, aes(
                x = sum(c(other_1, other_2)), xend = sum(c(sel_1, sel_2)),
                y = unif_feat_feat, yend = unif_feat_feat
            ), linetype = 2)
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = sum(c(other_1, other_2)), y = unif_feat_feat,
                    colour = "other cell expression"
                )
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = sum(c(sel_1, sel_2)), y = unif_feat_feat,
                    colour = "selected cell expression"
                )
            )
            pl <- pl + ggplot2::scale_colour_manual(
                name = "expression source", values = cols
            )
            pl <- pl + ggplot2::facet_wrap(~unif_int, scales = facet_scales)
            pl <- pl + ggplot2::labs(x = "interactions", y = "feat-feat")
        } else {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_bw()
            pl <- pl + ggplot2::geom_segment(data = subDT, aes(
                x = sum(c(other_1, other_2)), xend = sum(c(sel_1, sel_2)),
                y = unif_int, yend = unif_int
            ), linetype = 2)
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = sum(c(other_1, other_2)), y = unif_int,
                    colour = "other cell expression"
                )
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = sum(c(sel_1, sel_2)), y = unif_int,
                    colour = "selected cell expression"
                )
            )
            pl <- pl + ggplot2::scale_colour_manual(
                name = "expression source", values = cols
            )
            pl <- pl + ggplot2::facet_wrap(
                ~unif_feat_feat,
                scales = facet_scales
            )
            pl <- pl + ggplot2::labs(x = "feat-feat", y = "interactions")
        }
    }

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}







#' @title plotCombineICF
#' @name plotCombineICF
#' @description Create visualization for combined (pairwise) ICF scores
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param combIcfObject ICFscores, output from combineInteractionChangedFeats()
#' @param selected_interactions interactions to show
#' @param selected_feat_to_feat pairwise feature combinations to show
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or feats with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param colors vector with two colors to use
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' g_icf <- findInteractionChangedFeats(g,
#'     cluster_column = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17"), nr_permutations = 10
#' )
#'
#' combIcfObject <- combineInteractionChangedFeats(g_icf)
#'
#' plotCombineICF(
#'     gobject = g, combIcfObject = combIcfObject,
#'     selected_feat_to_feat = c("Btbd17--Ccnd2", "Btbd17--Gna12"),
#'     selected_interactions = "1--8"
#' )
#' @export
plotCombineICF <- function(
        gobject,
        combIcfObject,
        selected_interactions = NULL,
        selected_feat_to_feat = NULL,
        detail_plot = TRUE,
        simple_plot = FALSE,
        simple_plot_facet = c("interaction", "feats"),
        facet_scales = "fixed",
        facet_ncol = length(selected_feat_to_feat),
        facet_nrow = length(selected_interactions),
        colors = c("#9932CC", "#FF8C00"),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotCombineICF") {
    plotCombineInteractionChangedFeats(
        gobject = gobject,
        combIcfObject = combIcfObject,
        selected_interactions = selected_interactions,
        selected_feat_to_feat = selected_feat_to_feat,
        detail_plot = detail_plot,
        simple_plot = simple_plot,
        simple_plot_facet = simple_plot_facet,
        facet_scales = facet_scales,
        facet_ncol = facet_ncol,
        facet_nrow = facet_nrow,
        colors = colors,
        show_plot = show_plot,
        return_plot = return_plot,
        save_plot = save_plot,
        save_param = save_param,
        default_save_name = default_save_name
    )
}












# * ####
# cell communication plots ####

#' @title plotCombineCellCellCommunication
#' @name plotCombineCellCellCommunication
#' @description Create visualization for combined (pairwise) cell proximity
#' gene scores
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param combCCcom combined communcation scores, output from combCCcom()
#' @param selected_LR selected ligand-receptor pair
#' @param selected_cell_LR selected cell-cell interaction pair for
#' ligand-receptor pair
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or genes with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param colors vector with two colors to use
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' comScores <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = c("Gm19935", "2900040C04Rik", "Ccnd2"),
#'     feat_set_2 = c("9630013A20Rik", "Gna12", "Btbd17")
#' )
#'
#' exprCC <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935", feat_set_2 = "9630013A20Rik"
#' )
#' spatialCC <- spatCellCellcom(
#'     gobject = g, cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935", feat_set_2 = "9630013A20Rik", verbose = "a lot"
#' )
#'
#' combCCcom <- combCCcom(spatialCC = spatialCC, exprCC = exprCC)
#'
#' plotCombineCellCellCommunication(
#'     gobject = g, combCCcom = combCCcom,
#'     selected_LR = c("Gm19935-9630013A20Rik"), selected_cell_LR = c("1--1")
#' )
#' @export
plotCombineCellCellCommunication <- function(
        gobject,
        combCCcom,
        selected_LR = NULL,
        selected_cell_LR = NULL,
        detail_plot = TRUE,
        simple_plot = FALSE,
        simple_plot_facet = c("interaction", "genes"),
        facet_scales = "fixed",
        facet_ncol = length(selected_LR),
        facet_nrow = length(selected_cell_LR),
        colors = c("#9932CC", "#FF8C00"),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotCombineCellCellCommunication") {
    # data.table variables
    LR_comb <- LR_cell_comb <- lig_expr <- lig_expr_spat <- rec_expr <-
        rec_expr_spat <- LR_expr <- LR_expr_spat <- NULL

    ## check validity
    if (is.null(selected_cell_LR) | is.null(selected_LR)) {
        stop("You need to provide a selection of cell-cell interactions
            and genes-genes to plot")
    }

    subDT <- combCCcom[
        LR_comb %in% selected_LR & LR_cell_comb %in% selected_cell_LR
    ]

    # order interactions and gene-to-gene according to input
    subDT[, LR_comb := factor(LR_comb, levels = selected_LR)]
    subDT[, LR_cell_comb := factor(LR_cell_comb, levels = selected_cell_LR)]

    if (simple_plot == FALSE) {
        pl <- ggplot2::ggplot()
        pl <- pl + ggplot2::theme_bw()

        if (detail_plot == TRUE) {
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(x = 0, y = lig_expr, colour = "overall cell expression"),
                shape = 1
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = 0, y = lig_expr_spat,
                    colour = "spatial cell expression"
                ), shape = 1
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(x = rec_expr, y = 0, colour = "overall cell expression"),
                shape = 1
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = rec_expr_spat, y = 0,
                    colour = "spatial cell expression"
                ), shape = 1
            )
        }

        pl <- pl + ggplot2::geom_point(
            data = subDT,
            aes(x = rec_expr, y = lig_expr, colour = "overall cell expression"),
            size = 2
        )
        pl <- pl + ggplot2::geom_point(
            data = subDT,
            aes(
                x = rec_expr_spat, y = lig_expr_spat,
                colour = "spatial cell expression"
            ), size = 2
        )
        pl <- pl + ggplot2::geom_segment(data = subDT, aes(
            x = rec_expr, xend = rec_expr_spat,
            y = lig_expr, yend = lig_expr_spat
        ), linetype = 2)
        pl <- pl + ggplot2::labs(
            x = paste(subDT$receptor, subDT$rec_cell_type, sep = " in "),
            y = paste(subDT$ligand, subDT$lig_cell_type, sep = " in ")
        )
        pl <- pl + ggplot2::scale_colour_manual(
            name = "expression source", values = colors
        )
        pl <- pl + ggplot2::facet_wrap(~ LR_comb + LR_cell_comb,
            nrow = facet_nrow, ncol = facet_ncol,
            scales = facet_scales
        )
    } else {
        simple_plot_facet <- match.arg(
            arg = simple_plot_facet, choices = c("interaction", "genes")
        )

        if (simple_plot_facet == "interaction") {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_bw()
            pl <- pl + ggplot2::geom_segment(data = subDT, aes(
                x = LR_expr, xend = LR_expr_spat,
                y = LR_comb, yend = LR_comb
            ), linetype = 2)
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = LR_expr, y = LR_comb,
                    colour = "overall cell expression"
                )
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = LR_expr_spat, y = LR_comb,
                    colour = "spatial cell expression"
                )
            )
            pl <- pl + ggplot2::scale_colour_manual(
                name = "expression source", values = colors
            )
            pl <- pl + ggplot2::facet_wrap(~LR_cell_comb, scales = "fixed")
            pl <- pl + ggplot2::labs(x = "interactions", y = "gene-gene")
            pl
        } else {
            pl <- ggplot2::ggplot()
            pl <- pl + ggplot2::theme_bw()
            pl <- pl + ggplot2::geom_segment(data = subDT, aes(
                x = LR_expr, xend = LR_expr_spat,
                y = LR_cell_comb, yend = LR_cell_comb
            ), linetype = 2)
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = LR_expr, y = LR_cell_comb,
                    colour = "overall cell expression"
                )
            )
            pl <- pl + ggplot2::geom_point(
                data = subDT,
                aes(
                    x = LR_expr_spat, y = LR_cell_comb,
                    colour = "spatial cell expression"
                )
            )
            pl <- pl + ggplot2::scale_colour_manual(
                name = "expression source", values = colors
            )
            pl <- pl + ggplot2::facet_wrap(~LR_comb, scales = facet_scales)
            pl <- pl + ggplot2::labs(x = "gene-gene", y = "interactions")
        }
    }

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}



#' @title plotCombineCCcom
#' @name plotCombineCCcom
#' @description Create visualization for combined (pairwise) cell proximity
#' gene scores
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param combCCcom combined communcation scores, output from combCCcom()
#' @param selected_LR selected ligand-receptor pair
#' @param selected_cell_LR selected cell-cell interaction pair for
#' ligand-receptor pair
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or genes with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param colors vector with two colors to use
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' comScores <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = c("Gm19935", "2900040C04Rik", "Ccnd2"),
#'     feat_set_2 = c("9630013A20Rik", "Gna12", "Btbd17")
#' )
#'
#' exprCC <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935", feat_set_2 = "9630013A20Rik"
#' )
#' spatialCC <- spatCellCellcom(
#'     gobject = g, cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935", feat_set_2 = "9630013A20Rik", verbose = "a lot"
#' )
#'
#' combCCcom <- combCCcom(spatialCC = spatialCC, exprCC = exprCC)
#'
#' plotCombineCCcom(
#'     gobject = g, combCCcom = combCCcom,
#'     selected_LR = c("Gm19935-9630013A20Rik"), selected_cell_LR = c("1--1")
#' )
#' @export
plotCombineCCcom <- function(
        gobject,
        combCCcom,
        selected_LR = NULL,
        selected_cell_LR = NULL,
        detail_plot = TRUE,
        simple_plot = FALSE,
        simple_plot_facet = c("interaction", "genes"),
        facet_scales = "fixed",
        facet_ncol = length(selected_LR),
        facet_nrow = length(selected_cell_LR),
        colors = c("#9932CC", "#FF8C00"),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotCombineCCcom") {
    plotCombineCellCellCommunication(
        gobject = gobject,
        combCCcom = combCCcom,
        selected_LR = selected_LR,
        selected_cell_LR = selected_cell_LR,
        detail_plot = detail_plot,
        simple_plot = simple_plot,
        simple_plot_facet = simple_plot_facet,
        facet_scales = facet_scales,
        facet_ncol = facet_ncol,
        facet_nrow = facet_nrow,
        colors = colors,
        show_plot = show_plot,
        return_plot = return_plot,
        save_plot = save_plot,
        save_param = save_param,
        default_save_name = default_save_name
    )
}



#' @title plotCCcomHeatmap
#' @name plotCCcomHeatmap
#' @description Plots heatmap for ligand-receptor communication scores in
#' cell-cell interactions
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_params
#' @param comScores communinication scores from \code{\link{exprCellCellcom}}
#' or \code{\link{spatCellCellcom}}
#' @param selected_LR selected ligand-receptor combinations
#' @param selected_cell_LR selected cell-cell combinations for ligand-receptor
#' combinations
#' @param show_LR_names show ligand-receptor names
#' @param show_cell_LR_names show cell-cell names
#' @param show values to show on heatmap
#' @param cor_method correlation method used for clustering
#' @param aggl_method agglomeration method used by hclust
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' comScores <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = c("Gm19935", "2900040C04Rik", "Ccnd2"),
#'     feat_set_2 = c("9630013A20Rik", "Gna12", "Btbd17")
#' )
#'
#' plotCCcomHeatmap(gobject = g, comScores = comScores, show_plot = TRUE)
#' @export
plotCCcomHeatmap <- function(
        gobject,
        comScores,
        selected_LR = NULL,
        selected_cell_LR = NULL,
        show_LR_names = TRUE,
        show_cell_LR_names = TRUE,
        show = c("PI", "LR_expr", "log2fc"),
        cor_method = c("pearson", "kendall", "spearman"),
        aggl_method = c(
            "ward.D", "ward.D2", "single", "complete", "average",
            "mcquitty", "median", "centroid"
        ),
        gradient_color = NULL,
        gradient_style = c("divergent", "sequential"),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotCCcomHeatmap") {
    # get parameters
    cor_method <- match.arg(
        cor_method,
        choices = c("pearson", "kendall", "spearman")
    )
    aggl_method <- match.arg(aggl_method, choices = c(
        "ward.D", "ward.D2", "single", "complete",
        "average", "mcquitty", "median", "centroid"
    ))


    # data.table variables
    LR_comb <- LR_cell_comb <- NULL

    # plot method
    if (!is.null(selected_LR) & !is.null(selected_cell_LR)) {
        selDT <- comScores[LR_comb %in% selected_LR & LR_cell_comb %in%
            selected_cell_LR]
    } else if (!is.null(selected_LR)) {
        selDT <- comScores[LR_comb %in% selected_LR]
    } else if (!is.null(selected_cell_LR)) {
        selDT <- comScores[LR_cell_comb %in% selected_cell_LR]
    } else {
        selDT <- comScores
    }

    # creat matrix
    show <- match.arg(show, choices = c("PI", "LR_expr", "log2fc"))
    selDT_d <- data.table::dcast.data.table(
        selDT, LR_cell_comb ~ LR_comb,
        value.var = show, fill = 0
    )
    selDT_m <- dt_to_matrix(selDT_d)

    ## cells
    corclus_cells_dist <- stats::as.dist(
        1 - cor_flex(x = t_flex(selDT_m), method = cor_method)
    )
    hclusters_cells <- stats::hclust(
        d = corclus_cells_dist, method = aggl_method
    )
    clus_names <- rownames(selDT_m)
    names(clus_names) <- seq_along(clus_names)
    clus_sort_names <- clus_names[hclusters_cells$order]
    selDT[, LR_cell_comb := factor(LR_cell_comb, clus_sort_names)]

    ## genes
    corclus_genes_dist <- stats::as.dist(
        1 - cor_flex(x = selDT_m, method = cor_method)
    )
    hclusters_genes <- stats::hclust(
        d = corclus_genes_dist, method = aggl_method
    )
    clus_names <- colnames(selDT_m)
    names(clus_names) <- seq_along(clus_names)
    clus_sort_names <- clus_names[hclusters_genes$order]
    selDT[, LR_comb := factor(LR_comb, clus_sort_names)]



    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_raster(data = selDT, aes_string(
        x = "LR_cell_comb",
        y = "LR_comb", fill = show
    ))

    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank()
    )
    if (show_LR_names == TRUE) {
        pl <- pl + ggplot2::theme(
            axis.text.y = element_text(),
            axis.ticks.y = element_line()
        )
    }
    if (show_cell_LR_names == TRUE) {
        pl <- pl + ggplot2::theme(
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
            axis.ticks.x = element_line()
        )
    }
    pl <- pl + set_default_color_continuous_CCcom_heatmap(
        colors = gradient_color,
        instrs = instructions(gobject),
        style = gradient_style
    )
    pl <- pl + ggplot2::labs(x = "cell-cell", y = "ligand-receptor")

    return(plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}



#' @title plotCCcomDotplot
#' @name plotCCcomDotplot
#' @description Plots dotplot for ligand-receptor communication scores in
#' cell-cell interactions
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_params
#' @param comScores communication scores from \code{\link{exprCellCellcom}}
#' or \code{\link{spatCellCellcom}}
#' @param selected_LR selected ligand-receptor combinations
#' @param selected_cell_LR selected cell-cell combinations for ligand-receptor
#' combinations
#' @param show_LR_names show ligand-receptor names
#' @param show_cell_LR_names show cell-cell names
#' @param cluster_on values to use for clustering of cell-cell and
#' ligand-receptor pairs
#' @param cor_method correlation method used for clustering
#' @param aggl_method agglomeration method used by hclust
#' @param dot_color_gradient character. continuous colors to use. palette to
#' use or vector of colors to use (minimum of 2).
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' comScores <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = c("Gm19935", "2900040C04Rik", "Ccnd2"),
#'     feat_set_2 = c("9630013A20Rik", "Gna12", "Btbd17")
#' )
#'
#' plotCCcomDotplot(gobject = g, comScores = comScores, show_plot = TRUE)
#' @export
plotCCcomDotplot <- function(
        gobject,
        comScores,
        selected_LR = NULL,
        selected_cell_LR = NULL,
        show_LR_names = TRUE,
        show_cell_LR_names = TRUE,
        cluster_on = c("PI", "LR_expr", "log2fc"),
        cor_method = c("pearson", "kendall", "spearman"),
        aggl_method = c(
            "ward.D", "ward.D2", "single", "complete", "average",
            "mcquitty", "median", "centroid"
        ),
        dot_color_gradient = NULL,
        gradient_style = c("divergent", "sequential"),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotCCcomDotplot") {
    # get parameters
    cor_method <- match.arg(
        cor_method,
        choices = c("pearson", "kendall", "spearman")
    )
    aggl_method <- match.arg(aggl_method, choices = c(
        "ward.D", "ward.D2", "single", "complete",
        "average", "mcquitty", "median", "centroid"
    ))


    # data.table variables
    LR_comb <- LR_cell_comb <- sd <- NULL

    # plot method
    if (!is.null(selected_LR) & !is.null(selected_cell_LR)) {
        selDT <- comScores[LR_comb %in% selected_LR & LR_cell_comb %in%
            selected_cell_LR]
    } else if (!is.null(selected_LR)) {
        selDT <- comScores[LR_comb %in% selected_LR]
    } else if (!is.null(selected_cell_LR)) {
        selDT <- comScores[LR_cell_comb %in% selected_cell_LR]
    } else {
        selDT <- comScores
    }

    # creat matrix
    cluster_on <- match.arg(cluster_on, choices = c("PI", "LR_expr", "log2fc"))
    selDT_d <- data.table::dcast.data.table(
        selDT, LR_cell_comb ~ LR_comb,
        value.var = cluster_on, fill = 0
    )
    selDT_m <- dt_to_matrix(selDT_d)

    # remove zero variance
    sd_rows <- apply(selDT_m, 1, sd)
    sd_rows_zero <- names(sd_rows[sd_rows == 0])
    if (length(sd_rows_zero) > 0) {
        selDT_m <- selDT_m[
            !rownames(selDT_m) %in% sd_rows_zero,
        ]
    }

    sd_cols <- apply(selDT_m, 2, sd)
    sd_cols_zero <- names(sd_cols[sd_cols == 0])
    if (length(sd_cols_zero) > 0) {
        selDT_m <- selDT_m[
            , !colnames(selDT_m) %in% sd_cols_zero
        ]
    }



    ## cells
    corclus_cells_dist <- stats::as.dist(
        1 - cor_flex(x = t_flex(selDT_m), method = cor_method)
    )
    hclusters_cells <- stats::hclust(
        d = corclus_cells_dist, method = aggl_method
    )
    clus_names <- rownames(selDT_m)
    names(clus_names) <- seq_along(clus_names)
    clus_sort_names <- clus_names[hclusters_cells$order]
    selDT[, LR_cell_comb := factor(LR_cell_comb, clus_sort_names)]

    ## genes
    corclus_genes_dist <- stats::as.dist(
        1 - cor_flex(x = selDT_m, method = cor_method)
    )
    hclusters_genes <- stats::hclust(
        d = corclus_genes_dist, method = aggl_method
    )
    clus_names <- colnames(selDT_m)
    names(clus_names) <- seq_along(clus_names)
    clus_sort_names <- clus_names[hclusters_genes$order]
    selDT[, LR_comb := factor(LR_comb, clus_sort_names)]



    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_point(data = selDT, aes_string(
        x = "LR_cell_comb",
        y = "LR_comb", size = "pvalue", color = "log2fc"
    ))
    pl <- pl + ggplot2::theme_classic()
    if (show_LR_names == TRUE) {
        pl <- pl + ggplot2::theme(
            axis.text.y = element_text(),
            axis.ticks.y = element_line()
        )
    }
    if (show_cell_LR_names == TRUE) {
        pl <- pl + ggplot2::theme(
            axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
            axis.ticks.x = element_line()
        )
    }
    pl <- pl + set_default_color_continuous_CCcom_dotplot(
        colors = dot_color_gradient,
        instrs = instructions(gobject),
        midpoint = NULL,
        style = gradient_style,
        type = "fill"
    )
    pl <- pl + ggplot2::scale_size_continuous(range = c(5, 0.5)) +
        set_default_color_continuous_CCcom_dotplot(
            colors = dot_color_gradient,
            instrs = instructions(gobject),
            midpoint = NULL,
            style = gradient_style,
            type = "color"
        )
    pl <- pl + ggplot2::labs(x = "cell-cell", y = "ligand-receptor")

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}



#' @title plotRankSpatvsExpr
#' @name plotRankSpatvsExpr
#' @description Plots dotplot to compare ligand-receptor rankings from
#' spatial and expression information
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param combCC combined communication scores from \code{\link{combCCcom}}
#' @param expr_rnk_column column with expression rank information to use
#' @param spat_rnk_column column with spatial rank information to use
#' @param gradient_midpoint numeric. default = 10. midpoint of colors
#' @param dot_color_gradient character. continuous colors to use. palette to
#' use or vector of colors to use (minimum of 2).
#' @param gradient_style either 'divergent' (midpoint is used in color scaling)
#' or 'sequential' (scaled based on data range)
#' @param midpoint deprecated
#' @param size_range size ranges of dotplot
#' @param xlims x-limits, numerical vector of 2
#' @param ylims y-limits, numerical vector of 2
#' @param selected_ranks numerical vector, will be used to print out the
#' percentage of top spatial ranks are recovered
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' exprCC <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935", feat_set_2 = "9630013A20Rik"
#' )
#' spatialCC <- spatCellCellcom(
#'     gobject = g, cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935", feat_set_2 = "9630013A20Rik", verbose = "a lot",
#'     random_iter = 10
#' )
#'
#' combCC <- combCCcom(spatialCC = spatialCC, exprCC = exprCC)
#'
#' plotRankSpatvsExpr(gobject = g, combCC = combCC)
#' @export
plotRankSpatvsExpr <- function(
        gobject,
        combCC,
        expr_rnk_column = "LR_expr_rnk",
        spat_rnk_column = "LR_spat_rnk",
        dot_color_gradient = NULL,
        midpoint = deprecated(),
        gradient_midpoint = 10,
        gradient_style = c("divergent", "sequential"),
        size_range = c(0.01, 1.5),
        xlims = NULL,
        ylims = NULL,
        selected_ranks = c(1, 10, 20),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotRankSpatvsExpr") {
    # deprecate
    if (GiottoUtils::is_present(midpoint)) {
        deprecate_warn(
            "4.0.0",
            "Giotto::plotRankSpatvsExpr(midpoint = )",
            "Giotto::plotRankSpatvsExpr(gradient_midpoint = )"
        )
        gradient_midpoint <- midpoint
    }


    # data.table variables
    spt_rank <- variable <- value <- NULL

    total_rnks <- max(unique(combCC[[expr_rnk_column]]))

    rnk_list <- list()
    spt_list <- list()
    for (rnk in seq_len(total_rnks)) {
        mytab <- table(cut(
            sort(combCC[get(expr_rnk_column) == rnk][[
                spat_rnk_column
            ]]),
            breaks = seq(0, total_rnks, 1),
            labels = seq_len(total_rnks)
        ))
        rnk_list[[rnk]] <- mytab
        spt_list[[rnk]] <- names(mytab)
    }

    rnk_res <- data.table::as.data.table(do.call("rbind", rnk_list))
    rnk_res[, spt_rank := seq_len(total_rnks)]

    rnk_res_m <- data.table::melt.data.table(rnk_res, id.vars = "spt_rank")
    rnk_res_m[, spt_rank := as.numeric(spt_rank)]
    rnk_res_m[, variable := as.numeric(variable)]

    rnk_res_m[, diff := variable - spt_rank]

    for (i in selected_ranks) {
        perc_recovered <- 100 * (sum(rnk_res_m[abs(diff) < i]$value) /
            sum(rnk_res_m$value))
        cat(
            "for top ", i, " expression ranks, you recover ",
            round(perc_recovered, 2), "% of the highest spatial rank"
        )
    }


    # full plot
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(
        axis.text = element_blank()
    )
    pl <- pl + ggplot2::geom_point(
        data = rnk_res_m,
        ggplot2::aes(x = variable, y = spt_rank, size = value, color = value)
    )
    pl <- pl + set_default_color_continuous_CCcom_dotplot(
        colors = dot_color_gradient,
        instrs = instructions(gobject),
        midpoint = gradient_midpoint,
        style = gradient_style,
        type = "color",
        guide = guide_legend(title = "")
    )
    pl <- pl + ggplot2::scale_size_continuous(
        range = size_range, guide = "none"
    )
    pl <- pl + ggplot2::labs(x = "expression rank", y = "spatial rank")

    if (!is.null(xlims)) {
        pl <- pl + xlim(xlims)
    }

    if (!is.null(ylims)) {
        pl <- pl + ylim(ylims)
    }

    return(plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}




#' @title Create recovery plot
#' @name .plotRecovery_sub
#' @description Plots recovery plot to compare ligand-receptor rankings from
#' spatial and expression information
#' @param combCC combined communinication scores from \code{\link{combCCcom}}
#' @param first_col first column to use
#' @param second_col second column to use
#' @returns ggplot
#' @keywords internal
.plotRecovery_sub <- function(
        combCC,
        first_col = "LR_expr_rnk",
        second_col = "LR_spat_rnk") {
    # data.table variables
    concord <- perc <- not_concord <- secondrank <- secondrank_perc <- NULL

    mergeDT_filt <- combCC[get(first_col) == 1]

    mymat <- matrix(data = NA, nrow = max(combCC[[second_col]]), ncol = 2)
    for (i in seq_len(max(combCC[[second_col]]))) {
        mergeDT_filt[, concord := ifelse(get(second_col) <= i, "yes", "no")]
        mytable <- table(mergeDT_filt$concord)

        matching <- mytable["yes"]
        if (is.na(matching)) matching <- 0
        mymat[i, 1] <- matching

        non_matching <- mytable["no"]
        if (is.na(non_matching)) non_matching <- 0
        mymat[i, 2] <- non_matching
    }

    mymatDT <- data.table::as.data.table(mymat)
    colnames(mymatDT) <- c("concord", "not_concord")
    mymatDT[, perc := 100 * (concord / (concord + not_concord))]
    mymatDT[, secondrank := seq_len(nrow(mymatDT))]
    mymatDT[, secondrank_perc := (secondrank / max(secondrank)) * 100]

    # percentage explained
    perc_explained <- mymatDT[, sum(perc)] / (100 * nrow(mymat))
    cat("percentage explained = ", perc_explained)


    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()
    pl <- pl + ggplot2::geom_point(
        data = mymatDT,
        aes(x = secondrank_perc, y = perc)
    )
    pl <- pl + ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 100))
    pl <- pl + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 100))
    pl <- pl + ggplot2::geom_abline(slope = 1, intercept = 0, color = "blue")

    return(pl)
}





#' @title plotRecovery
#' @name plotRecovery
#' @description Plots recovery plot to compare ligand-receptor rankings from
#' spatial and expression information
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param combCC combined communication scores from \code{\link{combCCcom}}
#' @param expr_rnk_column column with expression rank information to use
#' @param spat_rnk_column column with spatial rank information to use
#' @param ground_truth what to consider as ground truth (default: spatial)
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' exprCC <- exprCellCellcom(g,
#'     cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935", feat_set_2 = "9630013A20Rik"
#' )
#' spatialCC <- spatCellCellcom(
#'     gobject = g, cluster_column = "leiden_clus",
#'     feat_set_1 = "Gm19935", feat_set_2 = "9630013A20Rik", verbose = "a lot",
#'     random_iter = 10
#' )
#'
#' combCC <- combCCcom(spatialCC = spatialCC, exprCC = exprCC)
#'
#' plotRecovery(gobject = g, combCC = combCC)
#' @export
plotRecovery <- function(
        gobject,
        combCC,
        expr_rnk_column = "exprPI_rnk",
        spat_rnk_column = "spatPI_rnk",
        ground_truth = c("spatial", "expression"),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotRecovery") {
    ground_truth <- match.arg(
        ground_truth,
        choices = c("spatial", "expression")
    )


    if (ground_truth == "spatial") {
        pl <- .plotRecovery_sub(
            combCC = combCC,
            first_col = spat_rnk_column,
            second_col = expr_rnk_column
        )
        pl <- pl + ggplot2::labs(
            x = "% expression rank included",
            y = "% highest spatial rank recovered"
        )
    } else if (ground_truth == "expression") {
        pl <- .plotRecovery_sub(
            combCC = combCC,
            first_col = expr_rnk_column,
            second_col = spat_rnk_column
        )
        pl <- pl + ggplot2::labs(
            x = "% spatial rank included",
            y = "% highest expression rank recovered"
        )
    }

    return(plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}









# * ####
# cell proximity spatplots ####

#' @title cellProximitySpatPlot2D
#' @name cellProximitySpatPlot2D
#' @description Visualize 2D cell-cell interactions according to spatial
#' coordinates in ggplot mode
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_spatnet_params
#' @inheritParams plot_params
#' @param spat_loc_name spatial locations to use
#' @param interaction_name cell-cell interaction name
#' @param cluster_column cluster column with cell clusters
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param show_other_cells decide if show cells not in network
#' @param show_network show spatial network of selected cells
#' @param show_other_network show spatial network of not selected cells
#' @param network_color color of spatial network
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param point_size_select size of selected points
#' @param point_select_border_col border color of selected points
#' @param point_select_border_stroke stroke size of selected points
#' @param point_size_other size of other points
#' @param point_alpha_other opacity of other points
#' @param point_other_border_col border color of other points
#' @param point_other_border_stroke stroke size of other points
#' @returns ggplot
#' @details Description of parameters.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' g <- createSpatialGrid(g, sdimx_stepsize = 5, sdimy_stepsize = 5)
#' x <- cellProximityEnrichment(g, cluster_column = "leiden_clus")
#'
#' cellProximitySpatPlot2D(
#'     gobject = g, cluster_column = "leiden_clus",
#'     interaction_name = x
#' )
#' @export
cellProximitySpatPlot2D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        spat_loc_name = NULL,
        interaction_name = NULL,
        cluster_column = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        cell_color = NULL,
        cell_color_code = NULL,
        color_as_factor = TRUE,
        show_other_cells = FALSE,
        show_network = FALSE,
        show_other_network = FALSE,
        network_color = NULL,
        spatial_network_name = "Delaunay_network",
        show_grid = FALSE,
        grid_color = NULL,
        spatial_grid_name = "spatial_grid",
        coord_fix_ratio = 1,
        show_legend = TRUE,
        point_size_select = 2,
        point_select_border_col = "black",
        point_select_border_stroke = 0.05,
        point_size_other = 1,
        point_alpha_other = 0.3,
        point_other_border_col = "lightgrey",
        point_other_border_stroke = 0.01,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "cellProximitySpatPlot2D") {
    if (is.null(interaction_name)) {
        stop("you need to specific at least one interaction name, run
             cellProximityEnrichment")
    }


    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # get information from all slots
    cell_locations <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table",
        copy_obj = FALSE
    )

    spatial_grid <- getSpatialGrid(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        name = spatial_grid_name,
        return_grid_Obj = FALSE
    )

    cell_metadata <- getCellMetadata(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        output = "data.table",
        copy_obj = TRUE
    )

    spatial_network <- annotateSpatialNetwork(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        spatial_network_name = spatial_network_name,
        cluster_column = cluster_column
    )


    # data.table variables
    unified_int <- sdimx_begin <- sdimy_begin <- sdimx_end <- sdimy_end <-
        x_start <- x_end <- y_start <- y_end <- cell_ID <- NULL

    cell_IDs_to_keep <- unique(c(
        spatial_network[unified_int %in% interaction_name]$to,
        spatial_network[unified_int %in% interaction_name]$from
    ))

    if (show_other_cells) {
        CellType <- strsplit(interaction_name, "--")
        all_cell_IDs <- cell_metadata[cell_metadata[[
            cluster_column
        ]] == CellType[[1]][1] |
            cell_metadata[[cluster_column]] == CellType[[1]][2], ]$cell_ID
        other_cell_IDs <- setdiff(all_cell_IDs, cell_IDs_to_keep)
    }


    # annotated cell data ()
    if (nrow(cell_metadata) == 0) {
        cell_locations_metadata <- cell_locations
    } else {
        cell_locations_metadata <- merge(
            cell_locations, cell_metadata,
            by = "cell_ID"
        )
    }


    # first 2 dimensions need to be defined
    if (is.null(sdimx) | is.null(sdimy)) {
        message("first and second dimension need to be defined, default is
                first 2")
        sdimx <- "sdimx"
        sdimy <- "sdimy"
    }

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()

    if (!is.null(spatial_network) & show_network == TRUE) {
        if (is.null(network_color)) network_color <- "red"
        if (show_other_network) {
            pl <- pl + ggplot2::geom_segment(
                data = spatial_network[!unified_int %in% interaction_name],
                aes(
                    x = sdimx_begin, y = sdimy_begin, xend = sdimx_end,
                    yend = sdimy_end
                ),
                color = "lightgrey", size = 0.5, alpha = 0.5
            )
        }
        pl <- pl + ggplot2::geom_segment(
            data = spatial_network[unified_int %in% interaction_name],
            aes(
                x = sdimx_begin, y = sdimy_begin, xend = sdimx_end,
                yend = sdimy_end
            ),
            color = network_color, size = 0.5, alpha = 0.5
        )
    }

    if (!is.null(spatial_grid) & show_grid == TRUE) {
        if (is.null(grid_color)) grid_color <- "black"
        pl <- pl + ggplot2::geom_rect(
            data = spatial_grid,
            aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end),
            color = grid_color, fill = NA
        )
    }

    # cell color default
    if (is.null(cell_color)) {
        cell_color <- "lightblue"
        pl <- pl + ggplot2::geom_point(
            data = cell_locations[!cell_ID %in% cell_IDs_to_keep],
            aes_string(x = sdimx, y = sdimy),
            show.legend = show_legend, shape = 21, fill = "lightgrey",
            size = point_size_other
        )
        pl <- pl + ggplot2::geom_point(
            data = cell_locations[cell_ID %in% cell_IDs_to_keep],
            aes_string(x = sdimx, y = sdimy),
            show.legend = show_legend, shape = 21, fill = cell_color,
            size = point_size_select
        )
        if (show_other_cells) {
            pl <- pl + ggplot2::geom_point(
                data = cell_locations[cell_ID %in% other_cell_IDs],
                aes_string(x = sdimx, y = sdimy),
                show.legend = show_legend, shape = 21, fill = cell_color,
                alpha = point_alpha_other,
                size = point_size_select * 0.5
            )
        }
    } else if (is.character(cell_color)) {
        if (cell_color %in% colnames(cell_locations_metadata)) {
            if (color_as_factor == TRUE) {
                factor_data <- factor(cell_locations_metadata[[cell_color]])
                cell_locations_metadata[[cell_color]] <- factor_data
            }

            pl <- pl + ggplot2::geom_point(
                data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
                aes_string(x = sdimx, y = sdimy),
                fill = "lightgrey", shape = 21, size = point_size_other,
                color = point_other_border_col,
                stroke = point_other_border_stroke
            )
            pl <- pl + ggplot2::geom_point(
                data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                aes_string(x = sdimx, y = sdimy, fill = cell_color),
                show.legend = show_legend, shape = 21, size = point_size_select,
                color = point_select_border_col,
                stroke = point_select_border_stroke
            )
            if (show_other_cells) {
                pl <- pl + ggplot2::geom_point(
                    data = cell_locations_metadata[cell_ID %in% other_cell_IDs],
                    aes_string(x = sdimx, y = sdimy, fill = cell_color),
                    show.legend = show_legend, shape = 21,
                    alpha = point_alpha_other,
                    size = point_size_select * 0.5
                )
            }



            if (!is.null(cell_color_code)) {
                pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
            } else if (color_as_factor == TRUE) {
                number_colors <- length(unique(factor_data))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instructions(gobject)
                )(n = number_colors)
                names(cell_color_code) <- unique(factor_data)
                pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
            } else if (color_as_factor == FALSE) {
                pl <- pl + set_default_color_continuous_cell(
                    colors = cell_color_code,
                    instrs = instructions(gobject),
                    midpoint = NULL,
                    style = "sequential",
                    type_default = list(
                        pal = c("blue", "red")
                    )
                )
            }
        } else {
            pl <- pl + ggplot2::geom_point(
                data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
                aes_string(x = sdimx, y = sdimy),
                show.legend = show_legend, shape = 21, fill = "lightgrey",
                size = point_size_other,
                color = point_other_border_col,
                stroke = point_other_border_stroke
            )
            pl <- pl + ggplot2::geom_point(
                data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                aes_string(x = sdimx, y = sdimy),
                show.legend = show_legend, shape = 21, fill = cell_color,
                size = point_size_select,
                color = point_select_border_col,
                stroke = point_select_border_stroke
            )
        }
    }

    pl <- pl + ggplot2::theme_bw() + ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)
    )

    # fix coord ratio
    if (!is.null(coord_fix_ratio)) {
        pl <- pl + ggplot2::coord_fixed(ratio = coord_fix_ratio)
    }

    pl <- pl + ggplot2::labs(x = "x coordinates", y = "y coordinates")

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}


#' @title cellProximitySpatPlot
#' @name cellProximitySpatPlot
#' @description Visualize 2D cell-cell interactions according to spatial
#' coordinates in ggplot mode
#' @param gobject giotto object
#' @inheritDotParams cellProximitySpatPlot2D -gobject
#' @returns ggplot
#' @details Description of parameters.
#' @export
#' @seealso  \code{\link{cellProximitySpatPlot2D}} and
#' \code{\link{cellProximitySpatPlot3D}} for 3D
cellProximitySpatPlot <- function(gobject, ...) {
    cellProximitySpatPlot2D(gobject = gobject, ...)
}


#' @title cellProximitySpatPlot3D
#' @name cellProximitySpatPlot3D
#' @description Visualize 3D cell-cell interactions according to spatial
#' coordinates in plotly mode
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_spatnet_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_params
#' @param interaction_name cell-cell interaction name
#' @param cluster_column cluster column with cell clusters
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimz')
#' @param show_other_cells decide if show cells not in network
#' @param show_network show spatial network of selected cells
#' @param show_other_network show spatial network of not selected cells
#' @param network_color color of spatial network
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param point_size_select size of selected points
#' @param point_size_other size of other points
#' @param point_alpha_other opacity of other points
#' @param axis_scale scale of axis
#' @param x_ticks ticks on x-axis
#' @param y_ticks ticks on y-axis
#' @param z_ticks ticks on z-axis
#' @param custom_ratio custom ratio of axes
#' @param \dots additional parameters
#' @returns plotly
#' @details Description of parameters.
#' @export
cellProximitySpatPlot3D <- function(
        gobject,
        interaction_name = NULL,
        cluster_column = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        sdimz = "sdimz",
        cell_color = NULL,
        cell_color_code = NULL,
        color_as_factor = TRUE,
        show_other_cells = TRUE,
        show_network = TRUE,
        show_other_network = FALSE,
        network_color = NULL,
        spatial_network_name = "Delaunay_network",
        show_grid = FALSE,
        grid_color = NULL,
        spatial_grid_name = "spatial_grid",
        show_legend = TRUE,
        point_size_select = 4,
        point_size_other = 2,
        point_alpha_other = 0.5,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        z_ticks = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "cellProximitySpatPlot3D",
        ...) {
    if (is.null(sdimz)) {
        pl <- .cellProximityVisPlot_2D_plotly(
            gobject = gobject,
            interaction_name = interaction_name,
            cluster_column = cluster_column,
            sdimx = sdimx,
            sdimy = sdimy,
            cell_color = cell_color,
            cell_color_code = cell_color_code,
            color_as_factor = color_as_factor,
            show_other_cells = show_other_cells,
            show_network = show_network,
            show_other_network = show_other_network,
            network_color = network_color,
            spatial_network_name = spatial_network_name,
            show_grid = show_grid,
            grid_color = grid_color,
            spatial_grid_name = spatial_grid_name,
            show_legend = show_legend,
            point_size_select = point_size_select,
            point_size_other = point_size_other,
            point_alpha_other = point_alpha_other,
            axis_scale = axis_scale,
            custom_ratio = custom_ratio,
            x_ticks = x_ticks,
            y_ticks = y_ticks,
            ...
        )
    } else {
        pl <- .cellProximityVisPlot_3D_plotly(
            gobject = gobject,
            interaction_name = interaction_name,
            cluster_column = cluster_column,
            sdimx = sdimx,
            sdimy = sdimy,
            sdimz = sdimz,
            cell_color = cell_color,
            cell_color_code = cell_color_code,
            color_as_factor = color_as_factor,
            show_other_cells = show_other_cells,
            show_network = show_network,
            show_other_network = show_other_network,
            network_color = network_color,
            spatial_network_name = spatial_network_name,
            show_grid = show_grid,
            grid_color = grid_color,
            spatial_grid_name = spatial_grid_name,
            show_legend = show_legend,
            point_size_select = point_size_select,
            point_size_other = point_size_other,
            point_alpha_other = point_alpha_other,
            axis_scale = axis_scale,
            custom_ratio = custom_ratio,
            x_ticks = x_ticks,
            y_ticks = y_ticks,
            z_ticks = z_ticks,
            ...
        )
    }

    # output plot
    return(GiottoVisuals::plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}
