# ------------------------------------------------------------------------- #
# This script contains reusable dummy documentation / templates for
# commonly used params.
#
# Use the @inheritParams tag to use these templates in documentation
# ------------------------------------------------------------------------- #

# No @title to prevent .Rd generation
# No @noRd tags can be used on these dummy documentations, otherwise @inheritParams
# cannot be used

# Note that these dummy documentations WILL be flagged with warnings when building
# the documentation, but this should be fine.








# Data Access ####


#' data_access_params
#'
#' @name data_access_params
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param return_uniques return unique nesting names (ignores if final object
#' exists/is correct class)
#' @param output what format in which to get information (e.g. "data.table")
#' @param set_defaults set default spat_unit and feat_type. Change to FALSE
#' only when
#' expression and spat_info are not expected to exist.
#' @param copy_obj whether to deep copy/duplicate when getting the object
#' (default = TRUE)
#' @param initialize (default = FALSE) whether to initialize the gobject before
#' returning
#' @returns list
#' @keywords internal
NULL




# Read Functions ####

#' read_data_params
#' @name read_data_params
#' @param data_list (nested) list of input data to read
#' @param default_spat_unit (optional) default spat_unit to use
#' @param default_feat_type (optional) default feat_type to use
#' @param provenance (optional) provenance information
#' @param verbose be verbose
#' @returns list
#' @keywords internal
NULL




# Plotting ####



#' Params documentation template: plot_cell_params
#' @name plot_cell_params
#' @param cell_color character. what to color cells by (e.g. metadata col or
#' spatial enrichment col)
#' @param color_as_factor logical. convert color column to factor. discrete
#' colors are used when this is TRUE. continuous colors when FALSE.
#' @param cell_color_code character. discrete colors to use. palette to use or
#' named vector of colors
#' @param cell_color_gradient character. continuous colors to use. palette to
#' use or vector of colors to use (minimum of 2).
#' @returns plot
#' @keywords internal
NULL

#' Params documentation template: plot_feat_params
#' @name plot_feat_params
#' @param feats_color_code code to color the provided features
#' @param feat_shape_code code to shape the provided feature types
#' @returns plot
#' @keywords internal
NULL

#' Params documentation template: plot_poly_params
#' @name plot_poly_params
#' @param show_polygon overlay polygon information (e.g. cell shape)
#' @param use_overlap use polygon and feature coordinates overlap results
#' @param polygon_feat_type feature type associated with polygon information
#' @param polygon_color color for polygon border
#' @param polygon_bg_color color for polygon background (overruled by
#' polygon_fill)
#' @param polygon_fill character. what to color to fill polgyons by (e.g.
#' metadata col or spatial enrichment col)
#' @param polygon_fill_gradient polygon fill gradient colors given in order
#' from low to high
#' @param polygon_fill_gradient_midpoint value to set as gradient midpoint
#' (optional). If left as \code{NULL}, the median value detected will be chosen
#' @param polygon_fill_gradient_style either 'divergent' (midpoint is used in
#' color scaling) or 'sequential' (scaled based on data range)
#' @param polygon_fill_as_factor is fill color a factor
#' @param polygon_fill_code code to color the fill column
#' @param polygon_alpha alpha of polygon
#' @param polygon_line_size line width of the polygon's outline
#' @returns plot
#' @keywords internal
NULL

#' Params documentation template: plot_dimred_params
#' @name plot_dimred_params
#' @param dim_reduction_to_use character. dimension reduction to use
#' @param dim_reduction_name character. dimension reduction name
#' @param dim1_to_use numeric. dimension to use on x-axis
#' @param dim2_to_use numeric. dimension to use on y-axis
#' @param dim3_to_use numeric. dimension to use on z-axis
#' @param dim_point_shape point with border or not (border or no_border)
#' @param dim_point_size size of points in dim. reduction space
#' @param dim_point_alpha transparancy of point in dim. reduction space
#' @param dim_point_border_col border color of points in dim. reduction space
#' @param dim_point_border_stroke border stroke of points in dim. reduction
#' space
#' @returns plot
#' @keywords internal
NULL

#' Params documentation template: plot_nn_net_params
#' @name plot_nn_net_params
#' @param show_NN_network logical. Show underlying NN network
#' @param nn_network_to_use character. type of NN network to use (kNN vs sNN)
#' @param network_name character. name of NN network to use, if
#' show_NN_network = TRUE
#' @param nn_network_name character. name of NN network to use, if
#' show_NN_network = TRUE
#' @param network_color color of NN network
#' @param nn_network_alpha column to use for alpha of the edges
#' @returns plot
#' @keywords internal
NULL

#' Params documentation template: plot_spatnet_params
#' @name plot_spatnet_params
#' @param show_spatial_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spat_network_name name of spatial network to use
#' @param spat_network_color color of spatial network
#' @param spatial_network_color color of spatial network
#' @param spat_network_alpha alpha of spatial network
#' @returns plot
#' @keywords internal
NULL

#' Params documentation template: plot_spatenr_params
#' @name plot_spatenr_params
#' @param spat_enr_names character. names of spatial enrichment results to
#' include
#' @returns plot
#' @keywords internal
NULL

#' Params documentation template: plot_image_params
#' @name plot_image_params
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image or multiple images with group_by
#' @param largeImage_name name of a giottoLargeImage or multiple images with
#' group_by
#' @returns plot
#' @keywords internal
NULL


#' Params documentation template: plot_params
#'
#' @name plot_params
#'
#' @param group_by character. Create multiple plots based on cell annotation
#' column
#' @param group_by_subset character. subset the group_by factor column
#'
#' @param gradient_midpoint numeric. midpoint for color gradient
#' @param gradient_style either 'divergent' (midpoint is used in color scaling)
#' or 'sequential' (scaled based on data range)
#' @param gradient_limits numeric vector with lower and upper limits
#' @param gradient_color character. continuous colors to use. palette to
#' use or vector of colors to use (minimum of 2).
#'
#' @param select_cell_groups select subset of cells/clusters based on
#' cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#'
#' @param show_other_cells display not selected cells
#' @param other_cell_color color for not selected cells
#' @param other_cell_alpha (0 to 1) alpha for not selected cells
#' @param other_cells_alpha (0 to 1) alpha for not selected cells
#' @param other_point_size point size for not selected cells
#'
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#'
#' @param center_point_size size of center points
#' @param center_point_border_col border color of center points
#' @param center_point_border_stroke border stroke size of center points
#'
#' @param label_size  size of labels
#' @param label_fontface font of labels
#'
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_alpha transparency of points
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#'
#' @param title character. title for plot, defaults to cell_color parameter
#' @param show_legend logical. show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @returns plot
#' @keywords internal
NULL


#' Params documentation template: plot_cow_params
#' @name plot_cow_params
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative heights of rows (e.g. c(1,2))
#' @param cow_rel_w cowplot param: relative widths of columns (e.g. c(1,2))
#' @param cow_align cowplot param: how to align
#' @returns cowplot
#' @keywords internal
NULL

#' Params documentation template: plot_output_params
#' @name plot_output_params
#' @param show_plot logical. show plot
#' @param return_plot logical. return ggplot object
#' @param save_plot logical. save the plot
#' @param save_param list of saving parameters, see
#' \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change,
#' change save_name in save_param
#' @returns plot
#' @keywords internal
NULL
